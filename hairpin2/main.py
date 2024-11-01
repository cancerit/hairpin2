# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import pysam
from hairpin2 import ref2seq as r2s, constants as c, helpers as h
import hairpin2
from statistics import mean, median, stdev
import argparse
import logging
import json
from itertools import tee
from typing import Literal
from collections.abc import Iterable


# N.B.
# pysam guards against:
# quality and seq length mismatch
# reference id is none
def flag_read_broad(
    read: pysam.AlignedSegment,
    vcf_start: int,
    min_mapqual: int,
    min_clipqual: int,
) -> int:
    invalid_flag = c.ValidatorFlags.CLEAR.value  # 0 - evaluates false

    try:
        mate_cig = str(read.get_tag('MC'))
    except KeyError:
        mate_cig = None
    if any(x is None for x in
            [read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig]):
        invalid_flag |= c.ValidatorFlags.READ_FIELDS_MISSING.value
    else:
        if not (read.flag & 0x2) or read.flag & 0xE00:
            invalid_flag |= c.ValidatorFlags.FLAG.value

        if read.mapping_quality < min_mapqual:
            invalid_flag |= c.ValidatorFlags.MAPQUAL.value

        if ('S' in read.cigarstring and  # type: ignore - program ensures can't be none
                mean(read.query_alignment_qualities) < min_clipqual):  # type: ignore - legit type issue here with pysam but I can't fix it
            invalid_flag |= c.ValidatorFlags.CLIPQUAL.value

        if (not (invalid_flag & c.ValidatorFlags.FLAG.value)
                and not (read.flag & 0x40)):
            read_range = range(read.reference_start,
                               read.reference_end)  # type: ignore - can't be none
            mate_range = range(read.next_reference_start,
                               r2s.ref_end_via_cigar(mate_cig,  # type: ignore
                                                     read.next_reference_start))
            ref_overlap = set(read_range).intersection(mate_range)
            if vcf_start in ref_overlap:
                invalid_flag |= c.ValidatorFlags.OVERLAP.value

    return invalid_flag


def flag_read_alt(
    read: pysam.AlignedSegment,
    vcf_start: int,
    vcf_stop: int,
    alt: str,
    mut_type: Literal['S', 'D', 'I'],
    min_basequal: int
) -> int:
    if mut_type not in ['S', 'D', 'I']:
        raise ValueError(
            'unsupported mut_type: {} - supports \'S\' (SUB) \'D\' (DEL) \'I\' (INS)'.format(mut_type))

    invalid_flag = c.ValidatorFlags.CLEAR.value

    if mut_type in ['S', 'I']:
        try:
            mut_pos = r2s.ref2querypos(read, vcf_start)
        except IndexError:
            invalid_flag |= c.ValidatorFlags.NOT_ALIGNED.value
        else:
            if mut_type == 'S':  # SUB
                # type: ignore - can't be none
                if read.query_sequence[mut_pos:mut_pos + len(alt)] != alt:
                    invalid_flag |= c.ValidatorFlags.NOT_ALT.value
                if any([bq < min_basequal
                        for bq
                        in read.query_qualities[mut_pos:mut_pos + len(alt)]]):  # type: ignore - can't be none
                    invalid_flag |= c.ValidatorFlags.BASEQUAL.value
            if mut_type == 'I':  # INS - mut_pos is position immediately before insertion
                if mut_pos + len(alt) > read.query_length:
                    invalid_flag |= c.ValidatorFlags.SHORT.value
                else:
                    mut_alns = [(q, r)
                                for q, r
                                in read.get_aligned_pairs()
                                if q in range(mut_pos + 1, mut_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        invalid_flag |= c.ValidatorFlags.BAD_OP.value
                    # type: ignore - can't be none
                    if read.query_sequence[mut_pos + 1:mut_pos + len(alt) + 1] != alt:
                        invalid_flag |= c.ValidatorFlags.NOT_ALT.value
    # DEL
    if mut_type == 'D':
        rng = list(range(vcf_start - 1, vcf_stop + 1))
        mut_alns = [q
                    for q, r
                    in read.get_aligned_pairs()
                    if r in rng]
        if len(mut_alns) != len(rng):
            invalid_flag |= c.ValidatorFlags.SHORT.value
        if (any([x is not None for x in mut_alns[1:-1]]) or
            any([x is None for x in [mut_alns[0], mut_alns[-1]]])):
                invalid_flag |= c.ValidatorFlags.BAD_OP.value

    return invalid_flag


# detect PCR duplicates previously missed due to slippage
# this implementation assumes that sorting on first element of each sublist
# is appropriate, per Peter's initial implementation.
# is an all against all comparison between all read lists more appropriate?
# and between pairs of readlists, why is comparing sorted pairs most appropriate?
def get_hidden_PCRdup_indices(
    readpair_ends: list[list[int]],
    max_span: int
) -> list[int]:
    dup_idcs: list[int] = []
    read_ends_sorted: list[tuple[int, list[int]]] = sorted([(i, sorted(l))
                                                            for i, l
                                                            in enumerate(readpair_ends)],
                                                           key=lambda x: x[1])
    # smallest first element. What was Peter's intention here?
    base_read_ends_list: list[list[int]] = [read_ends_sorted[0][1]]
    for i in range(1, len(read_ends_sorted)):
        comparison_read_ends = read_ends_sorted[i]
        max_diffs = []
        for sublist in base_read_ends_list:
            max_diffs.append(max([abs(x - y)
                                  for x, y
                                  in zip(sublist, comparison_read_ends[1])]))
        if all([x <= max_span for x in max_diffs]):
            # dups
            base_read_ends_list.append(comparison_read_ends[1])
            dup_idcs.append(comparison_read_ends[0])
        else:
            # read at i is not dup of reads in base_read_ends_list
            # start again, test read at i
            # against reads subsequent from i in ends_sorted
            base_read_ends_list = [comparison_read_ends[1]]
    return dup_idcs


def alt_filter_reads(
    vstart: int,
    vstop: int,
    alt: str,
    mut_type: str,
    region_reads_by_sample: dict[str, Iterable[pysam.AlignedSegment]],
    max_span: int = 6,
    min_basequal: int = 25
) -> list[pysam.AlignedSegment]:
    rrbs_filt: dict[str, list[pysam.AlignedSegment]] = {key: []
                                                        for key
                                                        in region_reads_by_sample}
    filtered_reads: list[pysam.AlignedSegment] = []

    for mut_sample, read_iter in region_reads_by_sample.items():
        sample_readpair_ends: list[list[int]] = []
        for read in read_iter:
            if not flag_read_alt(read,
                                 vstart,
                                 vstop,
                                 alt,
                                 mut_type,  # type: ignore - type checkers annoying about literals
                                 min_basequal):
                rrbs_filt[mut_sample].append(read)
                next_ref_end = r2s.ref_end_via_cigar(
                    str(read.get_tag('MC')),
                    read.next_reference_start)
                sample_readpair_ends.append([read.reference_start,
                                             read.reference_end,  # type: ignore - won't be unbound within program
                                             read.next_reference_start,
                                             next_ref_end])
        if len(rrbs_filt[mut_sample]) > 1:
            drop_idcs = get_hidden_PCRdup_indices(sample_readpair_ends,
                                                  max_span=max_span)
            filtered_reads = filtered_reads + [j
                                               for i, j
                                               in enumerate(rrbs_filt[mut_sample])
                                               if i not in drop_idcs]
    return filtered_reads


def is_variant_AL(
    mut_reads: Iterable[pysam.AlignedSegment],
    al_thresh: float = 0.93
) -> c.ALFilter:
    al_filt = c.ALFilter()
    aln_scores: list[float] = []

    for read in mut_reads:
        try:
            aln_scores.append(int(read.get_tag('AS')) / read.query_length)
        except KeyError:
            pass
    if len(aln_scores) != 0:
        al_filt.avg_as = median(aln_scores)
        al_filt.code = c.FiltCodes.ON_THRESHOLD.value
        if al_filt.avg_as <= al_thresh:
            al_filt.set()
    else:
        al_filt.code = c.FiltCodes.INSUFFICIENT_READS.value

    return al_filt


# per Peter's implementation
# can set hairpin for mutations nowhere near alignment start
# expose more ellis conditions as parameters?
def is_variant_AD(
    vstart: int,
    mut_reads: Iterable[pysam.AlignedSegment],
    edge_definition: float = 0.15,
    edge_clustering_threshold: float = 0.9,
    min_MAD_one_strand: int = 0,  # exclusive (and subsequent)
    min_sd_one_strand: float = 4,
    min_MAD_both_strand_weak: int = 2,
    min_sd_both_strand_weak: float = 2,
    min_MAD_both_strand_strong: int = 1,
    min_sd_both_strand_strong: float = 10,
    min_reads: int = 1  # inclusive
) -> c.ADFilter:

    ad_filt = c.ADFilter()
    # *l*engths of *a*lignment starts *to* *m*utant query positions
    la2ms_f: list[int] = []
    la2ms_r: list[int] = []
    near_start_f: list[bool] = []
    near_start_r: list[bool] = []

    for read in mut_reads:
        mut_qpos = r2s.ref2querypos(read, vstart)
        if read.flag & 0x10:
            # +1 to include last base in length
            la2m = read.query_alignment_end - mut_qpos + 1
            near_start_r.append(((la2m / read.query_alignment_length)
                                 <= edge_definition))
            la2ms_r.append(la2m)
        else:
            la2m = mut_qpos - read.query_alignment_start + 1
            near_start_f.append(((la2m / read.query_alignment_length)
                                 <= edge_definition))
            la2ms_f.append(la2m)

    # hairpin conditions from Ellis et al. 2020, Nature Protocols
    # sometimes reported as 2021
    if len(la2ms_f) <= min_reads and len(la2ms_r) <= min_reads:
        ad_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
    else:
        if len(la2ms_f) > min_reads:
            mad_f = max(la2ms_f) - min(la2ms_f)
            sd_f = stdev(la2ms_f)
            if len(la2ms_r) <= min_reads:
                if (((sum(near_start_f) / len(near_start_f)) < edge_clustering_threshold) and
                    mad_f > min_MAD_one_strand and
                        sd_f > min_sd_one_strand):
                    ad_filt.code = c.FiltCodes.SIXTYAI.value  # 60A(i)
                else:
                    ad_filt.code = c.FiltCodes.SIXTYAI.value
                    ad_filt.set()
        if len(la2ms_r) > min_reads:
            mad_r = max(la2ms_r) - min(la2ms_r)
            sd_r = stdev(la2ms_r)
            if len(la2ms_f) <= min_reads:
                if (((sum(near_start_r) / len(near_start_r)) < edge_clustering_threshold) and
                    mad_r > min_MAD_one_strand and
                        sd_r > min_sd_one_strand):
                    ad_filt.code = c.FiltCodes.SIXTYAI.value
                else:
                    ad_filt.code = c.FiltCodes.SIXTYAI.value
                    ad_filt.set()
        if len(la2ms_f) > min_reads and len(la2ms_r) > min_reads:
            frac_lt_thresh = (sum(near_start_f + near_start_r)
                              / (len(near_start_f) + len(near_start_r)))
            if (frac_lt_thresh < edge_clustering_threshold or
                # type: ignore
                (mad_f > min_MAD_both_strand_weak and mad_r > min_MAD_both_strand_weak and sd_f > min_sd_both_strand_weak and sd_r > min_sd_both_strand_weak) or
                # type: ignore
                (mad_f > min_MAD_both_strand_strong and sd_f > min_sd_both_strand_strong) or
                    (mad_r > min_MAD_both_strand_strong and sd_r > min_sd_both_strand_strong)):  # type: ignore
                ad_filt.code = c.FiltCodes.SIXTYBI.value  # 60B(i)
            else:
                ad_filt.code = c.FiltCodes.SIXTYBI.value
                ad_filt.set()

    return ad_filt


def test_record_all_alts(
    alignments: dict[str, pysam.AlignmentFile],
    vcf_rec: pysam.VariantRecord,
    min_mapqual: int,
    min_clipqual: int,
    min_basequal: int,
    max_span: int,
    al_thresh: float,
    edge_def: float,
    edge_frac: float,
    mos: int,
    sos: float,
    mbsw: int,
    sbsw: float,
    mbss: int,
    sbss: float,
    min_reads: int
) -> dict[str, c.Filters]:

    if vcf_rec.alts is None:
        raise c.NoAlts
    samples_w_mutants = [name
                         for name
                         in vcf_rec.samples
                         if vcf_rec.samples[name]["GT"] != (0, 0)]
    if len(samples_w_mutants) == 0:
        raise c.NoMutants

    region_reads_by_sample: dict[str, Iterable[pysam.AlignedSegment]] = {}
    for k, v in alignments.items():
        if k in samples_w_mutants:
            read_iter, test_iter = tee(v.fetch(vcf_rec.chrom,
                                               vcf_rec.start,
                                               (vcf_rec.start + 1)))
            try:
                next(test_iter)
            except StopIteration:
                continue
            else:
                broad_filtered_iter = (read
                                       for read
                                       in read_iter
                                       if not flag_read_broad(read,
                                                              vcf_rec.start,
                                                              min_mapqual,
                                                              min_clipqual))
                # doesn't check for overwrite
                region_reads_by_sample[k] = broad_filtered_iter

    filt_d = {}
    for alt in vcf_rec.alts:
        if (vcf_rec.rlen == len(alt)
                and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):
            mut_type = 'S'
        elif len(alt) < vcf_rec.rlen or alt == '.':  # DEL - DOES NOT SUPPORT <DEL> TYPE IDS
            mut_type = 'D'
        elif (vcf_rec.rlen == 1
                and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):  # INS - DOES NOT SUPPORT <INS> TYPE IDS
            mut_type = 'I'
        else:
            logging.warning('could not infer mutation type, POS={} REF={} ALT={}, skipping variant'.format(
                vcf_rec.pos, vcf_rec.ref, alt))
            continue
        alt_filt_reads: list = alt_filter_reads(vcf_rec.start,
                                                vcf_rec.stop,
                                                alt,
                                                mut_type,
                                                region_reads_by_sample,
                                                max_span,
                                                min_basequal)
        if len(alt_filt_reads) == 0:
            filt_d[alt] = c.Filters(c.ALFilter(code=c.FiltCodes.INSUFFICIENT_READS.value),
                                    c.ADFilter(code=c.FiltCodes.INSUFFICIENT_READS.value))
        else:
            filt_d[alt] = c.Filters(is_variant_AL(alt_filt_reads,
                                                  al_thresh),
                                    is_variant_AD(vcf_rec.start,
                                                  alt_filt_reads,
                                                  edge_def,
                                                  edge_frac,
                                                  mos,
                                                  sos,
                                                  mbsw,
                                                  sbsw,
                                                  mbss,
                                                  sbss,
                                                  min_reads))
    return filt_d


def main_cli() -> None:
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s',
                        datefmt='%I:%M:%S')

    parser = argparse.ArgumentParser(prog="hairpin2",
                                     description='cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6). See README for further explanation of parameters.')
    parser._optionals.title = 'info'
    parser.add_argument('-v',
                        '--version',
                        help='print version',
                        action='version',
                        version=hairpin2.__version__)
    req = parser.add_argument_group('mandatory')
    req.add_argument('-i',
                     '--vcf-in',
                     help="path to input VCF",
                     required=True)
    req.add_argument('-o',
                     '--vcf-out',
                     help="path to write output VCF",
                     required=True)
    req.add_argument('-a',
                     '--alignments',
                     help="list of paths to (S/B/CR)AMs (indicated by --format) for samples in input VCF, whitespace separated - (s/b/cr)ai expected in same directories",
                     nargs='+',
                     required=True)
    req.add_argument('-f',
                     "--format",
                     help="format of alignment files; s indicates SAM, b indicates BAM, and c indicates CRAM",
                     choices=["s", "b", "c"],
                     type=str,
                     required=True)
    opt_rv = parser.add_argument_group('read validation')
    opt_rv.add_argument('-mc',
                        '--min-clip-quality',
                        help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35, range: 0-93, exclusive',
                        type=int)
    opt_rv.add_argument('-mq',
                        '--min-mapping-quality',
                        help='discard reads with mapping quality below this value - default: 11, range: 0-60, exclusive',
                        type=int)
    opt_rv.add_argument('-mb',
                        '--min-base-quality',
                        help='discard reads with base quality at variant position below this value - default: 25, range: 0-93, exclusive',
                        type=int)
    opt_rv.add_argument('-ms',
                        '--max-read-span',
                        help='maximum +- position to use when detecting PCR duplicates. -1 will disable duplicate detection - default: 6, range: -1-, inclusive',
                        type=int)
    opt_fc = parser.add_argument_group('filter conditions')
    opt_fc.add_argument('-al',
                        '--al-filter-threshold',
                        help='ALF; threshold for median of read alignment score per base of all relevant reads, at and below which a variant is flagged as ALF - default: 0.93, range: 0-, inclusive',
                        type=float)
    opt_fc.add_argument('-ed',
                        '--edge-definition',
                        help='ADF; percentage of a read that is considered to be "the edge" for the purposes of assessing variant location distribution - default: 0.15, range: 0-0.99, inclusive',
                        type=float)
    opt_fc.add_argument('-ef',
                        '--edge-fraction',
                        help='ADF; percentage of variants must occur within EDGE_FRACTION of read edges to allow ADF flag - default: 0.15, range: 0-0.99, exclusive',
                        type=float)
    opt_fc.add_argument('-mos',
                        '--min-MAD-one-strand',
                        help='ADF; min range of distances between variant position and read start for valid reads when only one strand has sufficient valid reads for testing - default: 0, range: 0-, exclusive',
                        type=int)
    opt_fc.add_argument('-sos',
                        '--min-sd-one-strand',
                        help='ADF; min stdev of variant position and read start for valid reads when only one strand has sufficient valid reads for testing - default: 4, range: 0-, exclusive',
                        type=float)
    opt_fc.add_argument('-mbsw',
                        '--min-MAD-both-strand-weak',
                        help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbsw is true - default: 2, range: 0-, exclusive',
                        type=int)
    opt_fc.add_argument('-sbsw',
                        '--min-sd-both-strand-weak',
                        help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbsw is true- default: 2, range: 0-, exclusive',
                        type=float)
    opt_fc.add_argument('-mbss',
                        '--min-MAD-both-strand-strong',
                        help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbss is true - default: 1, range: 0-, exclusive',
                        type=int)
    opt_fc.add_argument('-sbss',
                        '--min-sd-both-strand-strong',
                        help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbss is true - default: 10, range: 0-, exclusive',
                        type=float)
    opt_fc.add_argument('-mr',
                        '--min-reads',
                        help='ADF; number of reads at and below which the hairpin filtering logic considers a strand to have insufficient reads for testing - default: 1, range: 0-',
                        type=int)
    proc = parser.add_argument_group('procedural')
    proc.add_argument('-r',
                      '--cram-reference',
                      help="path to FASTA format CRAM reference, overrides $REF_PATH and UR tags - ignored if --format is not CRAM")
    proc.add_argument('-m',
                      '--name-mapping',
                      help='map VCF sample names to alignment SM tags; useful if they differ',
                      metavar='VCF:aln',
                      nargs='+')
    proc.add_argument('-ji',
                      '--input-json',
                      help='path to JSON of input parameters, from which extended arguments will be loaded - overridden by arguments provided on command line',
                      type=str)
    proc.add_argument('-jo',
                      '--output-json',
                      help='log input arguments to JSON',
                      type=str)

    args = parser.parse_args()

    json_config: dict | None = None
    if args.input_json:
        logging.info(
            'args JSON provided, extended arguments will be loaded from JSON if not present on command line')
        try:
            with open(args.input_json, 'r') as f:
                json_config = json.load(f)
        except Exception as e:
            h.cleanup(msg='failed to open input JSON, reporting: {}'.format(e))

    # set arg defaults
    for k in vars(args).keys():
        if not vars(args)[k]:
            if (json_config and k
                in json_config.keys()
                    and k in c.DEFAULTS.keys()):
                setattr(args, k, json_config[k])
            elif k in c.DEFAULTS.keys():
                setattr(args, k, c.DEFAULTS[k])

    # test args are sensible, exit if not
    if not any([(0 <= args.min_clip_quality <= 93),
                (0 <= args.min_mapping_quality <= 60),
                (0 <= args.min_base_quality <= 93),
                (args.max_read_span >= -1),
                (args.al_filter_threshold >= 0),
                (0 <= args.edge_definition <= 0.99),
                (0 <= args.edge_fraction <= 0.99),
                (args.min_MAD_one_strand >= 0),
                (args.min_sd_one_strand >= 0),
                (args.min_MAD_both_strand_weak >= 0),
                (args.min_sd_both_strand_weak >= 0),
                (args.min_MAD_both_strand_strong >= 0),
                (args.min_sd_both_strand_strong >= 0),
                (args.min_reads >= 0)]):
        h.cleanup(msg='extended arg out of range, check helptext for ranges')

    try:
        vcf_in_handle = pysam.VariantFile(args.vcf_in)
    except Exception as e:
        h.cleanup(msg='failed to open VCF input, reporting: {}'.format(e))
    sample_names = list(vcf_in_handle.header.samples)  # type:ignore
    if len(set(sample_names)) != len(sample_names):
        h.cleanup(msg='duplicate sample names in VCF')
    sample_names: set[str] = set(sample_names)

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}
    match args.format:
        case "s":
            mode = "r"
            logging.info("SAM format specified")
        case "b":
            mode = "rb"
            logging.info("BAM format specified")
        case "c":
            mode = "rc"
            logging.info("CRAM format specified")
    for path in args.alignments:
        try:
            alignment = pysam.AlignmentFile(path,
                                            mode,  # type: ignore - argparse ensures not unbound
                                            reference_filename=(args.cram_reference
                                                                if args.cram_reference
                                                                and args.format == "c"
                                                                else None))
        except Exception as e:
            h.cleanup(
                msg='failed to read alignment file at {}, reporting: {}'.format(path, e))
        # grab the sample name from first SM field
        # in header field RG
        # type: ignore - program ensures not unbound
        alignment_sample_name = alignment.header.to_dict()['RG'][0]['SM']
        # type: ignore - program ensures not unbound
        vcf_sample_to_alignment_map[alignment_sample_name] = alignment
    if args.name_mapping:
        if len(args.name_mapping) > len(args.alignments):
            h.cleanup(msg="more name mappings than alignments provided")
        vcf_map_names = []
        alignment_map_names = []
        for pair in args.name_mapping:
            kv_split = pair.split(':')  # VCF:aln
            if len(kv_split) != 2:
                h.cleanup(
                    msg='name mapping misformatted, more than two elements in map string {}'.format(pair))
            vcf_map_names.append(kv_split[0])
            alignment_map_names.append(kv_split[1])
        if h.has_duplicates(vcf_map_names):
            h.cleanup(
                msg='duplicate VCF sample names provided to name mapping flag')
        if not set(vcf_map_names) <= sample_names:
            h.cleanup(
                msg="VCF sample names provided to name mapping flag are not equal to, or a subset of, VCF sample names as retrieved from VCF")
        if h.has_duplicates(alignment_map_names):
            h.cleanup(
                msg='duplicate aligment sample names provided to name mapping flag')
        if h.lists_not_equal(alignment_map_names,
                             vcf_sample_to_alignment_map.keys()):  # type: ignore - dicts are stable
            h.cleanup(
                msg='alignment sample names provided to name mapping flag do not match alignment SM tags')
        vcf_sample_to_alignment_map = {vcf_map_names[alignment_map_names.index(k)]: v
                                       for k, v
                                       in vcf_sample_to_alignment_map.items()}
    else:
        if not vcf_sample_to_alignment_map.keys() <= sample_names:
            h.cleanup(msg='alignment SM tags do not match VCF sample names: {}'.format(
                vcf_sample_to_alignment_map.keys() - sample_names))
    if sample_names != vcf_sample_to_alignment_map.keys():
        logging.info("alignments not provided for all VCF samples; {} will be ignored".format(
            sample_names - vcf_sample_to_alignment_map.keys()))

    # init output
    out_head = vcf_in_handle.header  # type:ignore
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}, using samples {}\">".format(
        args.al_filter_threshold, ', '.join(vcf_sample_to_alignment_map.keys())))
    out_head.add_line("##FILTER=<ID=ADF,Description=\"Variant arises from hairpin artefact, using samples {}\">".format(
        ', '.join(vcf_sample_to_alignment_map.keys())))
    out_head.add_line(
        "##INFO=<ID=ADF,Number=1,Type=String,Description=\"alt|code for each alt indicating hairpin filter decision code\">")
    out_head.add_line(
        "##INFO=<ID=ALF,Number=1,Type=String,Description=\"alt|code|score for each alt indicating AL filter conditions\">")

    try:
        vcf_out_handle = pysam.VariantFile(args.vcf_out, 'w', header=out_head)
    except Exception as e:
        h.cleanup(msg='failed to open VCF output, reporting: {}'.format(e))

    # write args once all verified
    if args.output_json:
        try:
            with open(args.output_json, "w") as output_json:
                json.dump(
                    {
                        k: vars(args)[k]
                        for k
                        in (vars(args).keys() - {'input_json', 'output_json', 'format'})
                    },
                    output_json, indent="")
        except Exception as e:
            h.cleanup(msg='failed to write output JSON, reporting: {}'.format(e))

    for record in vcf_in_handle.fetch():  # type: ignore - program ensures not unbound
        try:
            filter_d: dict[str, c.Filters] = test_record_all_alts(
                vcf_sample_to_alignment_map,
                record,
                args.min_mapping_quality,
                args.min_clip_quality,
                args.min_base_quality,
                args.max_read_span,
                args.al_filter_threshold,
                args.edge_definition,
                args.edge_fraction,
                args.min_MAD_one_strand,
                args.min_sd_one_strand,
                args.min_MAD_both_strand_weak,
                args.min_sd_both_strand_weak,
                args.min_MAD_both_strand_strong,
                args.min_sd_both_strand_strong,
                args.min_reads
            )
        except c.NoAlts:
            logging.warning('{0: <7}:{1: >12} ¦ no alts for this record'.format(
                record.chrom, record.pos))
        except c.NoMutants:
            logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit record alts'.format(
                record.chrom, record.pos))
        else:
            for alt, filter_bundle in filter_d.items():
                for filter in filter_bundle:
                    if filter.flag:
                        record.filter.add(filter.name)
                    record.info.update({filter.name: '|'.join(  # type: ignore - unclear what pysam wants
                        [alt, str(int(filter.flag)), str(filter.code)] +
                        ([str(filter.avg_as)] if filter.name == 'ALF' else [])
                    )})
            try:
                vcf_out_handle.write(record)  # type:ignore
            except Exception as e:
                h.cleanup(msg='failed to write to vcf, reporting: {}'.format(e))
    h.cleanup(c.EXIT_SUCCESS)
