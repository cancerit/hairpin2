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
from functools import partial
from typing import Literal


def validate_read(
    read: pysam.AlignedSegment,
    vcf_start: int,
    vcf_stop: int,
    alt: str,
    mut_type: Literal['S', 'D', 'I'],
    min_mapqual: int,
    min_clipqual: int,
    min_basequal: int,
) -> int:
    sup_char_alt = ['A', 'C', 'G', 'T', 'N', '*', '.']
    if any([b not in sup_char_alt for b in alt]):
        raise ValueError('unsupported character in alt: {} - supports {}'.format(alt, ', '.join(sup_char_alt)))
    if mut_type not in ['S', 'D', 'I']:
        raise ValueError('unsupported mut_type: {} - supports \'S\' (SUB) \'D\' (DEL) \'I\' (INS)'.format(mut_type))

    read_flag = c.ValidatorFlags.CLEAR.value

    try:
        mate_cig = read.get_tag('MC')
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
        read_flag |= c.ValidatorFlags.READ_FIELDS_MISSING.value
    else:
        if not (read.flag & 0x2) or read.flag & 0xE00:
            read_flag |= c.ValidatorFlags.FLAG.value

        if read.mapping_quality < min_mapqual:
            read_flag |= c.ValidatorFlags.MAPQUAL.value

        if ('S' in read.cigarstring and  # type: ignore
                mean(read.query_alignment_qualities) < min_clipqual):  # type: ignore
            read_flag |= c.ValidatorFlags.CLIPQUAL.value

        if mut_type == 'S':  # SUB
            try:
                mut_pos, mut_op = r2s.ref2querypos(read, vcf_start)
            except IndexError:
                read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
            else:
                if read.query_sequence[mut_pos:mut_pos + len(alt)] != alt:  # type: ignore
                    read_flag |= c.ValidatorFlags.NOT_ALT.value
                if any([bq < min_basequal for bq in read.query_qualities[mut_pos:mut_pos + len(alt)]]):  # type: ignore
                    read_flag |= c.ValidatorFlags.BASEQUAL.value
        elif mut_type == 'D':  # DEL - doesn't check for matches before and after...
            # this could error if read doesn't cover region (as could all)
            mut_alns = [q for q, r in read.get_aligned_pairs() if r in range(vcf_start, vcf_stop)]
            if any([x is not None for x in mut_alns]):
                read_flag |= c.ValidatorFlags.BAD_OP.value
        elif mut_type == 'I':  # INS
            try:
                prior_pos, _ = r2s.ref2querypos(read, vcf_start)
            except IndexError:
                read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
            else:
                if prior_pos + len(alt) > read.query_length:
                    read_flag |= c.ValidatorFlags.SHORT.value
                else:
                    mut_alns = [(q, r) for q, r in read.get_aligned_pairs() if q in range(prior_pos + 1, prior_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        read_flag |= c.ValidatorFlags.BAD_OP.value
                    if read.query_sequence[prior_pos + 1:prior_pos + len(alt) + 1] != alt:
                        read_flag |= c.ValidatorFlags.NOT_ALT.value

        if read_flag == c.ValidatorFlags.CLEAR.value and not (read.flag & 0x40):
            read_range = range(read.reference_start,
                               read.reference_end)
            mate_range = range(read.next_reference_start,
                               r2s.ref_end_via_cigar(mate_cig,
                                                     read.next_reference_start)
                               )
            ref_overlap = set(read_range).intersection(mate_range)
            if vcf_start in ref_overlap:
                read_flag |= c.ValidatorFlags.OVERLAP.value

    return read_flag


def test_variant(
    vcf_rec: pysam.VariantRecord,
    mutant_alignments: dict[str, pysam.AlignmentFile],
    alt: str,
    mut_type: str,
    al_thresh: float,
    max_span: int,
    position_fraction_thresh: float,
    read_validator: c.FlagReturn,
) -> c.Filters:

    hp_filt = c.HPFilter()
    al_filt = c.ALFilter()

    mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in mutant_alignments}
    mut_reads_log: dict[str, list[tuple]] = {key: [] for key in mutant_alignments}
    mut_read_pos_f: list[int] = []
    mut_read_pos_r: list[int] = []
    mut_read_fracs_f: list[float] = []
    mut_read_fracs_r: list[float] = []
    aln_scores: list[float] = []

    for mut_sample, alignment in mutant_alignments.items():
        read_iter, test_iter = tee(alignment.fetch(vcf_rec.chrom,
                                                   vcf_rec.start,
                                                   (vcf_rec.start + 1)))
        try:
            next(test_iter)
        except StopIteration:
            continue
        sample_readpair_ends = []
        read = None
        for read in read_iter:  # type: ignore
            read_flag = c.ValidatorFlags.CLEAR.value
            read_flag = read_validator(read=read,
                                       alt=alt,
                                       vcf_start=vcf_rec.start,
                                       vcf_stop=vcf_rec.stop,
                                       mut_type=mut_type)

            if read_flag == c.ValidatorFlags.CLEAR.value:
                mut_reads[mut_sample].append(read)
                sample_readpair_ends.append(
                            [read.reference_start,
                                read.reference_end,
                                read.next_reference_start,
                                r2s.ref_end_via_cigar(
                                                read.get_tag('MC'),
                                                read.next_reference_start)])  # type: ignore
            mut_reads_log[mut_sample].append((read.query_name, read_flag))
            del (read)
        if len(mut_reads[mut_sample]) > 1:
            sample_readpair_ends_sorted: list[list[int]] = sorted(list(map(
                                                        sorted,
                                                        sample_readpair_ends)))
            curr_ends = [sample_readpair_ends_sorted[0]]
            drop_idx = []
            for i in range(1, len(sample_readpair_ends_sorted)):
                max_spans = map(lambda sublist:
                                max(
                                    [abs(x - y)
                                        for x, y
                                        in zip(sublist,
                                               sample_readpair_ends_sorted[i])
                                     ]
                                ),
                                curr_ends)
                if all([x <= max_span for x in max_spans]):
                    curr_ends.append(sample_readpair_ends_sorted[i])
                    drop_idx.append(i)
                else:
                    curr_ends = [sample_readpair_ends_sorted[i]]
            mut_reads[mut_sample] = [j
                                     for i, j
                                     in enumerate(mut_reads[mut_sample])
                                     if i not in drop_idx]
    if all([len(x) == 0 for x in mut_reads.values()]):
        al_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
        hp_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
    else:
        for read_list in mut_reads.values():
            for read in read_list:
                mut_pos, _ = r2s.ref2querypos(read, vcf_rec.start)
                if read.flag & 0x10:
                    # 1-based position where start, idx 1, is alignment end
                    read_idx_wrt_aln = read.query_alignment_end - mut_pos
                    mut_read_fracs_r.append(read_idx_wrt_aln
                                            / read.query_alignment_length)
                    mut_read_pos_r.append(read_idx_wrt_aln)
                else:
                    read_idx_wrt_aln = mut_pos - read.query_alignment_start + 1
                    mut_read_fracs_f.append(read_idx_wrt_aln
                                            / read.query_alignment_length)
                    mut_read_pos_f.append(read_idx_wrt_aln)
                try:
                    aln_scores.append(read.get_tag('AS') / read.query_length)  # type:ignore
                except KeyError:
                    pass
        if len(aln_scores) != 0:
            al_filt.avg_as = median(aln_scores)
            al_filt.code = c.FiltCodes.ON_THRESHOLD.value
            if al_filt.avg_as <= al_thresh:
                al_filt.set()
        else:
            al_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
        # hairpin conditions from Ellis et al.
        if len(mut_read_pos_f) > 1 and not len(mut_read_pos_r) > 1:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
            if (
                ((sum([x <= position_fraction_thresh
                      for x
                      in mut_read_fracs_f]) / len(mut_read_pos_f)) < 0.9) and
                mad_f > 0 and
                    sd_f > 4):
                hp_filt.code = c.FiltCodes.SIXTYAI.value  # 60A(i)
            else:
                hp_filt.code = c.FiltCodes.SIXTYAI.value
                hp_filt.set()
        elif len(mut_read_pos_r) > 1 and not len(mut_read_pos_f) > 1:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            if (
                ((sum([x <= position_fraction_thresh
                      for x
                      in mut_read_fracs_r]) / len(mut_read_pos_r)) < 0.9) and
                mad_r > 0 and
                    sd_r > 4):
                hp_filt.code = c.FiltCodes.SIXTYAI.value
            else:
                hp_filt.code = c.FiltCodes.SIXTYAI.value
                hp_filt.set()
        elif len(mut_read_pos_f) > 1 and len(mut_read_pos_r) > 1:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            frac_lt_thresh = (sum([x <= position_fraction_thresh
                                  for x
                                  in mut_read_fracs_f + mut_read_fracs_r]) /
                              (len(mut_read_pos_f) + len(mut_read_pos_r)))
            if (frac_lt_thresh < 0.9 or
               (mad_f > 2 and mad_r > 2 and sd_f > 2 and sd_r > 2) or
               (mad_f > 1 and sd_f > 10) or
               (mad_r > 1 and sd_r > 10)):
                hp_filt.code = c.FiltCodes.SIXTYBI.value  # 60B(i)
            else:
                hp_filt.code = c.FiltCodes.SIXTYBI.value
                hp_filt.set()
        else:
            hp_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
    return c.Filters(al_filt, hp_filt)


def test_record_per_alt(
    alignments: dict[str, pysam.AlignmentFile],
    vcf_rec: pysam.VariantRecord,
    variant_tester: c.FiltReturn,
) -> dict[str, c.Filters]:

    if vcf_rec.alts is None:
        raise c.NoAlts
    samples_w_mutants = [name
                         for name
                         in vcf_rec.samples
                         if vcf_rec.samples[name]["GT"] != (0, 0)]
    if len(samples_w_mutants) == 0:
        raise c.NoMutants

    alignments_w_mutants = {k: v
                            for k, v
                            in alignments.items()
                            if k in samples_w_mutants}
    filt_d = {}
    for alt in vcf_rec.alts:
        if vcf_rec.rlen == len(alt) and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*'])):
            mut_type = 'S'
        elif len(alt) < vcf_rec.rlen or alt == '.':  # DEL - DOES NOT SUPPORT <DEL> TYPE IDS
            mut_type = 'D'
        elif vcf_rec.rlen == 1 and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*'])):  # INS - DOES NOT SUPPORT <INS> TYPE IDS
            mut_type = 'I'
        else:
            ## ERROR
            logging.warning('could not type mutation POS={} REF={} ALT={}, skipping alt'.format(vcf_rec.pos, vcf_rec.ref, alt))
            continue
        filt_d[alt] = variant_tester(vcf_rec, alignments_w_mutants, alt, mut_type)
    return filt_d


def main_cli() -> None:
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s',
                        datefmt='%I:%M:%S')

    parser = argparse.ArgumentParser(prog="hairpin2",
                                     description='cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6)')
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
    opt = parser.add_argument_group('extended')
    opt.add_argument('-al',
                     '--al-filter-threshold',
                     help='threshold for median of read alignment score per base of all relevant reads, below which a variant is flagged as ALF - default: 0.93',
                     type=float)
    opt.add_argument('-mc',
                     '--min-clip-quality',
                     help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35',
                     type=int)
    opt.add_argument('-mq',
                     '--min-mapping-quality',
                     help='discard reads with mapping quality below this value - default: 11',
                     type=int)
    opt.add_argument('-mb',
                     '--min-base-quality',
                     help='discard reads with base quality at variant position below this value - default: 25',
                     type=int)
    opt.add_argument('-ms',
                     '--max-read-span',
                     help='maximum +- position to use when detecting PCR duplicates - default: 6',
                     type=int)
    opt.add_argument('-pf',
                     '--position-fraction',
                     help='>90%% of variant must occur within POSITION_FRACTION of read edges to allow HPF flag - default: 0.15',
                     type=float)
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
        logging.info('args JSON provided, extended arguments will be loaded from JSON if not present on command line')
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
    h.test_options(args)

    primed_validate_read = partial(validate_read,
                                   min_mapqual=args.min_mapping_quality,
                                   min_clipqual=args.min_clip_quality,
                                   min_basequal=args.min_base_quality)

    primed_variant_tester = partial(test_variant,
                                    al_thresh=args.al_filter_threshold,
                                    max_span=args.max_read_span,
                                    position_fraction_thresh=args.position_fraction,
                                    read_validator=primed_validate_read)

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
                                            mode,
                                            reference_filename=(args.cram_reference
                                                                if args.cram_reference
                                                                and args.format == "c"
                                                                else None))
        except Exception as e:
            h.cleanup(msg='failed to read alignment file at {}, reporting: {}'.format(path, e))
        # grab the sample name from first SM field
        # in header field RG
        alignment_sample_name = alignment.header.to_dict()['RG'][0]['SM']
        vcf_sample_to_alignment_map[alignment_sample_name] = alignment
    if args.name_mapping:
        if len(args.name_mapping) > len(args.alignments):
            h.cleanup(msg="more name mappings than alignments provided")
        vcf_map_names = []
        alignment_map_names = []
        for pair in args.name_mapping:
            kv_split = pair.split(':')  # VCF:aln
            if len(kv_split) != 2:
                h.cleanup(msg='name mapping misformatted, more than two elements in map string {}'.format(pair))
            vcf_map_names.append(kv_split[0])
            alignment_map_names.append(kv_split[1])
        if h.has_duplicates(vcf_map_names):
            h.cleanup(msg='duplicate VCF sample names provided to name mapping flag')
        if not set(vcf_map_names) <= sample_names:
            h.cleanup(msg="VCF sample names provided to name mapping flag are not equal to, or a subset of, VCF sample names as retrieved from VCF")
        if h.has_duplicates(alignment_map_names):
            h.cleanup(msg='duplicate aligment sample names provided to name mapping flag')
        if h.lists_not_equal(alignment_map_names,
                             vcf_sample_to_alignment_map.keys()):
            h.cleanup(msg='alignment sample names provided to name mapping flag do not match alignment SM tags')
        vcf_sample_to_alignment_map = {vcf_map_names[alignment_map_names.index(k)]: v
                                       for k, v
                                       in vcf_sample_to_alignment_map.items()}
    else:
        if not vcf_sample_to_alignment_map.keys() <= sample_names:
            h.cleanup(msg='alignment SM tags do not match VCF sample names: {}'.format(vcf_sample_to_alignment_map.keys() - sample_names))
    if sample_names != vcf_sample_to_alignment_map.keys():
        logging.info("alignments not provided for all VCF samples; {} will be ignored".format(sample_names - vcf_sample_to_alignment_map.keys()))

    # init output
    out_head = vcf_in_handle.header  # type:ignore
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}, using samples {}\">".format(args.al_filter_threshold, ', '.join(vcf_sample_to_alignment_map.keys())))
    out_head.add_line("##FILTER=<ID=HPF,Description=\"Variant arises from hairpin artefact, using samples {}\">".format(', '.join(vcf_sample_to_alignment_map.keys())))
    out_head.add_line("##INFO=<ID=HPF,Number=1,Type=String,Description=\"alt|code for each alt indicating hairpin filter decision code\">")
    out_head.add_line("##INFO=<ID=ALF,Number=1,Type=String,Description=\"alt|code|score for each alt indicating AL filter conditions\">")

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

    for record in vcf_in_handle.fetch():  # type:ignore
        # need to test pysam's vcf record validation
        # e.g. what if start is after end
        try:
            filter_d: dict[str, c.Filters] = test_record_per_alt(
                alignments=vcf_sample_to_alignment_map,
                vcf_rec=record,
                variant_tester=primed_variant_tester
            )
        except c.NoAlts:
            logging.warning('{0: <7}:{1: >12} ¦ no alts for this record'.format(record.chrom, record.pos))
        except c.NoMutants:
            logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit record alts'.format(record.chrom, record.pos))
        else:
            for alt, filter_bundle in filter_d.items():
                for filter in filter_bundle:
                    if filter.flag:
                        record.filter.add(filter.name)
                    record.info.update({filter.name: '|'.join(
                                                      [alt] +
                                                      [str(f)
                                                          if type(f)
                                                          is not float
                                                          else str(round(f, 3))
                                                          for f in filter
                                                       ][2:]
                                                    )})

            try:
                vcf_out_handle.write(record)  # type:ignore
            except Exception as e:
                h.cleanup(msg='failed to write to vcf, reporting: {}'.format(e))
    h.cleanup(c.EXIT_SUCCESS)
