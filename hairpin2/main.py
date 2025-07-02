# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# pyright: reportUnusedCallResult=false
# pyright: reportAny=false
# pyright: reportExplicitAny=false

import pysam
from hairpin2 import ref2seq as r2s, constants as cnst, helpers as hlp
from hairpin2.filters import ADF, ALF, DVF, AnyResultCollection
import hairpin2
from statistics import mean
import argparse
import logging
import json
from itertools import tee, chain
from typing import Literal, Any
from collections.abc import Iterable
import sys
from tqdm import tqdm


def qc_read_broad(
    read: pysam.AlignedSegment,
    vcf_start: int,
    min_mapqual: int,
    min_clipqual: int,
) -> cnst.ValidatorFlags:
    """
    When testing variants, test a read for various general - i.e. not specific to a
    particular variant - features identifying the read as a poor source of support
    for a variant. Used to disqualify reads for use in testing a variant.
    
    does not check for the following, as pysam already guards against:
        - quality and seq length mismatch
        - reference id is none
    """
    invalid_flag = cnst.ValidatorFlags.CLEAR  # 0 - evaluates false

    try:
        mate_cig = str(read.get_tag('MC'))
    except KeyError:
        mate_cig = None
    if any(flg is None for flg in
            [read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig]):
        invalid_flag |= cnst.ValidatorFlags.READ_FIELDS_MISSING
    else:
        if not (read.flag & 0x2) or read.flag & 0xE00:
            invalid_flag |= cnst.ValidatorFlags.FLAG

        if read.mapping_quality < min_mapqual:
            invalid_flag |= cnst.ValidatorFlags.MAPQUAL

        if ('S' in read.cigarstring and  # pyright: ignore[reportOperatorIssue]
                mean(read.query_alignment_qualities) < min_clipqual):  # pyright: ignore[reportUnknownMemberType, reportArgumentType]
            invalid_flag |= cnst.ValidatorFlags.CLIPQUAL

        # avoid analysing both read1 and mate if they both cover the variant
        if (not (invalid_flag & cnst.ValidatorFlags.FLAG)
                and not (read.flag & 0x40)):
            read_range = range(read.reference_start,
                               read.reference_end)  # pyright: ignore[reportArgumentType]
            mate_range = range(read.next_reference_start,
                               r2s.ref_end_via_cigar(mate_cig,  # pyright: ignore[reportArgumentType]
                                                     read.next_reference_start))
            ref_overlap = set(read_range).intersection(mate_range)
            if vcf_start in ref_overlap:
                invalid_flag |= cnst.ValidatorFlags.OVERLAP

    return invalid_flag


def qc_read_alt_specific(
    read: pysam.AlignedSegment,
    vcf_start: int,
    vcf_stop: int,
    alt: str,
    mut_type: Literal['S', 'D', 'I'],
    min_basequal: int
) -> cnst.ValidatorFlags:
    """
    When testing a variant, test a read for various features specific to the variant
    at hand which would identify that read as a poor source of support for that
    variant. Used to disqualify reads for use in testing a variant.
    """
    if mut_type not in ['S', 'D', 'I']:
        raise ValueError(
            'unsupported mut_type: {} - supports \'S\' (SUB) \'D\' (DEL) \'I\' (INS)'.format(mut_type))
    if read.query_sequence is None:
        raise ValueError(
            'read must have query sequence'
        )
    if read.query_qualities is None:
        raise ValueError(
            'read must have query qualities'
        )

    invalid_flag = cnst.ValidatorFlags.CLEAR

    if mut_type in ['S', 'I']:
        try:
            mut_pos = r2s.ref2querypos(read, vcf_start)
        except ValueError:
            invalid_flag |= cnst.ValidatorFlags.NOT_ALIGNED
        else:
            if mut_type == 'S':  # SUB
                if read.query_sequence[mut_pos:mut_pos + len(alt)] != alt:
                    invalid_flag |= cnst.ValidatorFlags.NOT_ALT
                if any([bq < min_basequal
                        for bq
                        in read.query_qualities[mut_pos:mut_pos + len(alt)]]):
                    invalid_flag |= cnst.ValidatorFlags.BASEQUAL
            if mut_type == 'I':  # INS - mut_pos is position immediately before insertion
                if mut_pos + len(alt) > read.query_length:
                    invalid_flag |= cnst.ValidatorFlags.SHORT
                else:
                    mut_alns = [(q, r)
                                for q, r
                                in read.get_aligned_pairs()
                                if q in range(mut_pos + 1, mut_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        invalid_flag |= cnst.ValidatorFlags.BAD_OP
                    if read.query_sequence[mut_pos + 1:mut_pos + len(alt) + 1] != alt:
                        invalid_flag |= cnst.ValidatorFlags.NOT_ALT
    # DEL
    if mut_type == 'D':
        rng = list(range(vcf_start, vcf_stop + 1))
        mut_alns = [q
                    for q, r
                    in read.get_aligned_pairs()
                    if r in rng]
        if len(mut_alns) != len(rng):
            invalid_flag |= cnst.ValidatorFlags.SHORT
        if (any([x is not None for x in mut_alns[1:-1]]) or
            any([x is None for x in [mut_alns[0], mut_alns[-1]]])):
                invalid_flag |= cnst.ValidatorFlags.BAD_OP

    return invalid_flag


def test_record_all_alts(
    alignments: dict[str, pysam.AlignmentFile],
    vcf_rec: pysam.VariantRecord,
    min_mapqual: int,
    min_clipqual: int,
    min_basequal: int,
    adf_params: ADF.Params,
    alf_params: ALF.Params,
    dvf_params: DVF.Params
) -> dict[str, AnyResultCollection]:

    if vcf_rec.alts is None:
        raise cnst.NoAlts
    samples_w_mutants = [name
                         for name
                         in vcf_rec.samples
                         if vcf_rec.samples[name]["GT"] != (0, 0)]
    if len(samples_w_mutants) == 0:
        raise cnst.NoMutants

    region_reads_by_sample: dict[str, Iterable[pysam.AlignedSegment]] = {}
    for k, v in alignments.items():
        if k in samples_w_mutants:
            read_iter, test_iter = tee(v.fetch(vcf_rec.chrom,
                                               vcf_rec.start,
                                               (vcf_rec.start + 1)))
            try:
                _ = next(test_iter)
            except StopIteration:
                continue
            else:
                broad_qc_region_reads = (
                    read for read in read_iter
                    if not qc_read_broad(
                        read,
                        vcf_rec.start,
                        min_mapqual,
                        min_clipqual
                    )
                )
                # doesn't check for overwrite
                region_reads_by_sample[k] = broad_qc_region_reads

    # the ability to handle complex mutations would be a potentially interesting future feature
    # for extending to more varied artifacts
    alt_filters: dict[str, AnyResultCollection] = {}
    for alt in vcf_rec.alts:
        if (vcf_rec.rlen == len(alt)
                and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):
            mut_type = 'S'
        elif (len(alt) < vcf_rec.rlen
                and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):  # DEL - DOES NOT SUPPORT <DEL> TYPE IDS OR .
            mut_type = 'D'
        elif (vcf_rec.rlen == 1
                and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):  # INS - DOES NOT SUPPORT <INS> TYPE IDS
            mut_type = 'I'
        else:
            logging.warning('could not infer mutation type, POS={} REF={} ALT={}, skipping variant'.format(
                vcf_rec.pos, vcf_rec.ref, alt))
            continue

        alt_qc_region_reads: dict[str, list[pysam.AlignedSegment]] = {
            key: [] for key in region_reads_by_sample
        }
        for mut_sample, read_iter in region_reads_by_sample.items():
            alt_qc_region_reads[mut_sample] = [
                read for read in read_iter
                if not qc_read_alt_specific(
                    read,
                    vcf_rec.start,
                    vcf_rec.stop,
                    alt,
                    mut_type,
                    min_basequal
                )
            ]

        # instantiate filters to test the QC'd reads
        ad = ADF.Filter(fixed_params=adf_params)
        al = ALF.Filter(fixed_params=alf_params)
        dv = DVF.Filter(fixed_params=dvf_params)

        # run the tests
        dup_qc_region_reads, dv_result = dv.test(alt, alt_qc_region_reads)
        testing_reads = list(chain.from_iterable(dup_qc_region_reads.values()))
        testing_reads, al_result = al.test(alt, testing_reads)
        testing_reads, ad_result = ad.test(alt, vcf_rec.start, testing_reads)
        alt_filters[alt] = (al_result, ad_result, dv_result)
    return alt_filters


def main_cli() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s',
        datefmt='%I:%M:%S'
    ) # TODO: optionally log to file
    parser = argparse.ArgumentParser(
        prog="hairpin2",
        description='cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6). See README for further explanation of parameters.'
    )
    parser._optionals.title = 'info'
    parser.add_argument(
        '-v',
        '--version',
        help='print version',
        action='version',
        version=hairpin2.__version__
    )
    req = parser.add_argument_group('mandatory')
    req.add_argument(
        '-i',
        '--vcf-in',
        help="path to input VCF",
        required=True
    )
    req.add_argument(
        '-o',
        '--vcf-out',
        help="path to write output VCF",
        required=True
    )
    req.add_argument(
        '-a',
        '--alignments',
        help="list of paths to (S/B/CR)AMs (indicated by --format) for samples in input VCF, whitespace separated - (s/b/cr)ai expected in same directories",
        nargs='+',
        required=True
    )
    req.add_argument(
        '-f',
        "--format",
        help="format of alignment files; s indicates SAM, b indicates BAM, and c indicates CRAM",
        choices=["s", "b", "c"],
        type=str,
        required=True
    )
    opt_rv = parser.add_argument_group('read validation')
    opt_rv.add_argument(
        '-mc',
        '--min-clip-quality',
        help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35, range: 0-93, exclusive',
        type=int
    )
    opt_rv.add_argument(
        '-mq',
        '--min-mapping-quality',
        help='discard reads with mapping quality below this value - default: 11, range: 0-60, exclusive',
        type=int
    )
    opt_rv.add_argument(
        '-mb',
        '--min-base-quality',
        help='discard reads with base quality at variant position below this value - default: 25, range: 0-93, exclusive',
        type=int
    )
    opt_rv.add_argument(
        '-dw',
        '--duplication-window-size',
        help='maximum window size, in number of bases to use when detecting PCR duplicates. -1 will disable duplicate detection - default: 6, range: -1-, inclusive',
        type=int
    )  # TODO: consider keeping -ms --max-read-span as a deprecated option for backwards compat
    opt_fc = parser.add_argument_group('filter conditions')
    opt_fc.add_argument(
        '-al',
        '--al-filter-threshold',
        help='ALF; threshold for median of read alignment score per base of all relevant reads, at and below which a variant is flagged as ALF - default: 0.93, range: 0-, inclusive',
        type=float
    )
    opt_fc.add_argument(
        '-ed',
        '--edge-definition',
        help='ADF; percentage of a read that is considered to be "the edge" for the purposes of assessing variant location distribution - default: 0.15, range: 0-0.99, inclusive',
        type=float
    )
    opt_fc.add_argument(
        '-ef',
        '--edge-fraction',
        help='ADF; percentage of variants must occur within EDGE_FRACTION of read edges to allow ADF flag - default: 0.15, range: 0-0.99, exclusive',
        type=float
    )
    opt_fc.add_argument(
        '-mos',
        '--min-MAD-one-strand',
        help='ADF; min range of distances between variant position and read start for valid reads when only one strand has sufficient valid reads for testing - default: 0, range: 0-, exclusive',
        type=int
    )
    opt_fc.add_argument(
        '-sos',
        '--min-sd-one-strand',
        help='ADF; min stdev of variant position and read start for valid reads when only one strand has sufficient valid reads for testing - default: 4, range: 0-, exclusive',
        type=float
    )
    opt_fc.add_argument(
        '-mbsw',
        '--min-MAD-both-strand-weak',
        help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbsw is true - default: 2, range: 0-, exclusive',
        type=int
    )
    opt_fc.add_argument(
        '-sbsw',
        '--min-sd-both-strand-weak',
        help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbsw is true- default: 2, range: 0-, exclusive',
        type=float
    )
    opt_fc.add_argument(
        '-mbss',
        '--min-MAD-both-strand-strong',
        help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbss is true - default: 1, range: 0-, exclusive',
        type=int
    )
    opt_fc.add_argument(
        '-sbss',
        '--min-sd-both-strand-strong',
        help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbss is true - default: 10, range: 0-, exclusive',
        type=float
    )
    opt_fc.add_argument(
        '-mr',
        '--min-reads',
        help='ADF; number of reads at and below which the hairpin filtering logic considers a strand to have insufficient reads for testing - default: 1, range: 0-, inclusive',
        type=int
    )
    proc = parser.add_argument_group('procedural')
    proc.add_argument(
        '-r',
        '--cram-reference',
        help="path to FASTA format CRAM reference, overrides $REF_PATH and UR tags - ignored if --format is not CRAM",
        type=str
    )
    proc.add_argument(
        '-m',
        '--name-mapping',
        help='key to map samples in a multisample VCF to alignment/s provided to -a. Uses VCF sample names per VCF header and alignment SM tags. With multiple alignments to -a, accepts a space separated list of sample:SM pairs. With a single alignment, also accepts a comma separated string of one or more possible sample-of-interest names like TUMOR,TUMOUR',
        nargs='+'
    )
    proc.add_argument(
        '-ji',
        '--input-json',
        help='path to JSON of command line parameters, from which arguments will be loaded - overridden by arguments provided at runtime',
        type=str
    )
    proc.add_argument(
        '-jo',
        '--output-json',
        help='log command line arguments to JSON',
        type=str
    )
    proc.add_argument(
        '--no-prog',
        help='disable progress bar',
        action='store_true'
    )
    args = parser.parse_args()

    json_config: dict[str, Any] | None = None
    if args.input_json:
        logging.info(
            'args JSON provided, non-path and non-mandatory arguments will be loaded from JSON if not present on command line')
        try:
            with open(args.input_json, 'r') as f:
                json_config = json.load(f)
        except Exception as er:
            logging.error(f'failed to open input JSON, reporting: {er}')
            sys.exit(cnst.EXIT_FAILURE)

    # set arg defaults
    for k in vars(args).keys():
        if not vars(args)[k]:
            if (json_config and k
                in json_config.keys()):
                setattr(args, k, json_config[k])
            elif k in cnst.DEFAULTS.keys():
                setattr(args, k, cnst.DEFAULTS[k])

    # prepare args for recording to header
    arg_d: dict[str, Any] = vars(args)
    rec_args = sorted(arg_d.keys() - {"vcf_in", "vcf_out", "alignments", "input_json", "output_json", "format"})

    # test args are sensible, exit if not
    if not any([(0 <= args.min_clip_quality <= 93),
                (0 <= args.min_mapping_quality <= 60),
                (0 <= args.min_base_quality <= 93),
                (args.duplication_window_size >= -1),
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
        logging.error('extended arg out of range, check helptext for ranges')
        sys.exit(cnst.EXIT_FAILURE)

    # set up filter params based on inputs TODO: move validation there, rather than above validation check
    al_params = ALF.Params(
        args.al_filter_threshold
    )
    ad_params = ADF.Params(
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
    dv_params = DVF.Params(
        args.duplication_window_size
    )

    try:
        vcf_in_handle = pysam.VariantFile(args.vcf_in)
    except Exception as er:
        logging.error(f'failed to open input VCF, reporting {er}')
        sys.exit(cnst.EXIT_FAILURE)
    vcf_names: list[str] = list(vcf_in_handle.header.samples)
    if len(set(vcf_names)) != len(vcf_names):
        logging.error('duplicate sample names in input VCF')
        sys.exit(cnst.EXIT_FAILURE)

    sm_to_aln_map: dict[str, pysam.AlignmentFile] = {}
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
        case _:
            raise ValueError(f'Unknown mode {args.format} found in args.format')
    for path in args.alignments:
        try:
            alignment = pysam.AlignmentFile(path,
                                            mode,
                                            reference_filename=(args.cram_reference
                                                                if args.cram_reference
                                                                and args.format == "c"
                                                                else None))
        except Exception as er:
            logging.error(f'failed to read alignment file {path!r}, reporting {er}')
            sys.exit(cnst.EXIT_FAILURE)
        # grab the sample name from first SM field
        # in header field RG
        aln_sm = alignment.header.to_dict()['RG'][0]['SM']
        sm_to_aln_map[aln_sm] = alignment

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}
    if args.name_mapping:
        vcf_mapflag: list[str] = []
        alignment_mapflag: list[str] = []
        if len(args.name_mapping) <= len(args.alignments) and all(m.count(':') == 1 for m in args.name_mapping) and not any("," in m for m in args.name_mapping):
            for pair in args.name_mapping:
                kv_split = pair.split(':')  # VCF:aln
                vcf_mapflag.append(kv_split[0])
                alignment_mapflag.append(kv_split[1])

            if hlp.has_duplicates(vcf_mapflag):
                logging.error('duplicate VCF sample names provided to --name-mapping flag')
                sys.exit(cnst.EXIT_FAILURE)

            if not set(vcf_mapflag) <= set(vcf_names):
                logging.error(f'VCF sample names provided to --name-mapping flag {vcf_mapflag!r} are not equal to or a subset of VCF samples from input VCF {vcf_names!r}')
                sys.exit(cnst.EXIT_FAILURE)

            if hlp.has_duplicates(alignment_mapflag):
                logging.error(msg='duplicate aligment sample names provided to name mapping flag')
                sys.exit(cnst.EXIT_FAILURE)

            if not hlp.collection_equal(alignment_mapflag, sm_to_aln_map.keys()):
                logging.error(msg='SM tags provided to name mapping flag {} are not equal to SM tag list from alignment files {} - flag recieved {}'.format(
                        alignment_mapflag,
                        set(sm_to_aln_map.keys()),
                        args.name_mapping))
                sys.exit(cnst.EXIT_FAILURE)

            vcf_sample_to_alignment_map = {vcf_mapflag[alignment_mapflag.index(k)]: v
                                            for k, v
                                            in sm_to_aln_map.items()}
        elif not any(":" in m for m in args.name_mapping) and len(args.alignments) == len(args.name_mapping) == 1:
            if "," in args.name_mapping[0]:
                possible_sample_names = args.name_mapping[0].split(',')
            else:
                possible_sample_names = args.name_mapping  # list of len 1
            matches = [n for n in possible_sample_names if n in set(vcf_names)]
            if len(matches) > 1:
                logging.error(msg='More than one of the VCF sample names provided to name mapping flag match any sample names in input VCF {} - flag recieved {}'.format(
                        vcf_names,
                        args.name_mapping))
                sys.exit(cnst.EXIT_FAILURE)
            elif not matches:
                logging.error(msg='None of the VCF sample names provided to name mapping flag match any sample names in input VCF {} - flag recieved {}'.format(
                        vcf_names,
                        args.name_mapping))
                sys.exit(cnst.EXIT_FAILURE)
            else:
                logging.info('matched alignment to sample {} in VCF'.format(matches[0]))
                vcf_sample_to_alignment_map[matches[0]] = alignment  # pyright: ignore[reportPossiblyUnboundVariable] | since length of alignments == 1, can just reuse this variable

        else:
            logging.error(msg='name mapping misformatted, see helptext for expectations - flag recieved: {}'.format(args.name_mapping))
            sys.exit(cnst.EXIT_FAILURE)
    else:  # no name mapping flag
        if not sm_to_aln_map.keys() <= set(vcf_names):
            logging.error(
                msg='alignment SM tags {} are not equal to or a subset of VCF sample names {}'.format(
                    set(sm_to_aln_map.keys()),
                    set(vcf_names)
                )
            )
            sys.exit(cnst.EXIT_FAILURE)
        vcf_sample_to_alignment_map = sm_to_aln_map
    ## end name mapping handling

    if set(vcf_names) != vcf_sample_to_alignment_map.keys():
        logging.info(
            "alignments not provided for all VCF samples; {} will be ignored".format(
                set(vcf_names) - vcf_sample_to_alignment_map.keys()
            )
        )

    # init output  TODO: put these in filter modules as funcs
    out_head = vcf_in_handle.header
    out_head.add_line(
        f'##FILTER=<ID=ALF,Description="Median alignment score of reads reporting variant less than {args.al_filter_threshold}, using samples {vcf_sample_to_alignment_map.keys()!r}">'
    )
    out_head.add_line(
        f'##FILTER=<ID=ADF,Description="Variant arises from hairpin artefact, using samples {vcf_sample_to_alignment_map.keys()}">'
    )
    out_head.add_line(
        '##FILTER=<ID=DVF,Description="PLACEHOLDER">'
    )  # TODO: placeholder!
    out_head.add_line(
        '##INFO=<ID=ADF,Number=1,Type=String,Description="alt|code for each alt indicating hairpin filter decision code">'
    )
    out_head.add_line(
        '##INFO=<ID=ALF,Number=1,Type=String,Description="alt|code|score for each alt indicating AL filter conditions">'
    )
    out_head.add_line(
        '##INFO=<ID=DVF,Number=1,Type=String,Description="alt|code|score for each alt indicating DV filter conditions">'
    )  # TODO: placeholder!
    out_head.add_line(
        f'##hairpin2_version={hairpin2.__version__}'
    )
    out_head.add_line(
        f'##hairpin2_params=[{json.dumps({k: arg_d[k] for k in rec_args})}]'
    )

    try:
        vcf_out_handle = pysam.VariantFile(args.vcf_out, 'w', header=out_head)
    except Exception as e:
        logging.error(msg='failed to open VCF output, reporting: {}'.format(e))
        sys.exit(cnst.EXIT_FAILURE)
    for record in vcf_in_handle.fetch() if args.no_prog else tqdm(vcf_in_handle.fetch(), total=None, unit='records'):
        try:
            filter_d: dict[str, AnyResultCollection] = test_record_all_alts(
                vcf_sample_to_alignment_map,
                record,
                args.min_mapping_quality,
                args.min_clip_quality,
                args.min_base_quality,
                adf_params=ad_params,
                alf_params=al_params,
                dvf_params=dv_params
            )
        except cnst.NoAlts:
            logging.warning('{0: <7}:{1: >12} ¦ no alts for this record'.format(
                record.chrom, record.pos))
            sys.exit(cnst.EXIT_FAILURE)
        except cnst.NoMutants:
            logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit alts associated with this record'.format(
                record.chrom, record.pos))  # should this be recorded in the VCF?
            sys.exit(cnst.EXIT_FAILURE)
        else:
            for _, filter_bundle in filter_d.items():
                for filter in filter_bundle:
                    if filter.flag:
                        record.filter.add(filter.name)
                    if finfo := filter.getinfo():
                        record.info.update({filter.name: finfo})  # pyright: ignore[reportArgumentType]
            try:
                vcf_out_handle.write(record)
            except Exception as e:
                logging.error(msg='failed to write to vcf, reporting: {}'.format(e))
                sys.exit(cnst.EXIT_FAILURE)

    # write args once all verified
    if args.output_json:
        try:
            with open(args.output_json, "w") as output_json:
                json.dump(
                    {
                        k: arg_d[k]
                        for k
                        in rec_args
                    },
                    output_json, indent="")
        except Exception as e:
            logging.error(msg='failed to write output JSON, reporting: {}'.format(e))
            sys.exit(cnst.EXIT_FAILURE)

    logging.info('hairpin complete')
    sys.exit(cnst.EXIT_SUCCESS)
