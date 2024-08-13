import pysam
from hairpin2 import ref2seq as r2s, constants as c, helpers as h
from statistics import mean, median, stdev
import argparse
import logging
import json
from itertools import tee
from functools import partial

# CIGAR best retrieved from CG:B,I tag - implement in future
def validate_read(
    vcf_record: pysam.VariantRecord,
    read: pysam.AlignedSegment,
    min_mapqual: int,
    min_clipqual: int,
    min_basequal: int,
    alt: str
) -> int:
    read_flag = c.ValidatorFlags.CLEAR.value

    if not (read.flag & 0x2) or read.flag & 0xE00:
        read_flag |= c.ValidatorFlags.FLAG.value

    if read.mapping_quality < min_mapqual:
        read_flag |= c.ValidatorFlags.MAPQUAL.value

    try:
        mate_cig = read.get_tag('MC')
    except KeyError:
        mate_cig = None
    if any(x is None for x in [read.reference_start,
                            read.reference_end,
                            read.query_sequence,
                            read.query_qualities,
                            read.query_alignment_qualities,
                            read.cigarstring,
                            read.cigartuples,
                            mate_cig]):
        read_flag |= c.ValidatorFlags.READ_FIELDS_MISSING.value
    else:
        if ('S' in read.cigarstring and  # type: ignore
            mean(read.query_alignment_qualities) < min_clipqual):  # type: ignore
            read_flag |= c.ValidatorFlags.CLIPQUAL.value
        # First, check for sub
        try:
            mut_pos, mut_op = r2s.ref2querypos(read, vcf_record.start) # VCF 1-INDEXED, BAM 0-INDEXED (vcf_record.start = 0-indexed mutation position)
        except IndexError:
            read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
        else:
            if vcf_record.rlen == len(alt) == 1:
                if (mut_op not in [c.Ops.MATCH.value, c.Ops.DIFF.value]):
                    read_flag |= c.ValidatorFlags.BAD_OP.value
                if read.query_sequence[mut_pos] != alt:  # type: ignore
                    read_flag |= c.ValidatorFlags.NOT_ALT.value
                if read.query_qualities[mut_pos] < min_basequal:  # type: ignore
                        read_flag |= c.ValidatorFlags.BASEQUAL.value
            # Second, check whether length of read can accommodate size of indel
            elif (mut_pos + vcf_record.rlen > read.query_length or
                  mut_pos + len(alt) > read.query_length):
                read_flag |= c.ValidatorFlags.SHORT.value
            else:
                if len(alt) == 1:  # DEL
                    try:
                        mut_rng = list(map(lambda x: r2s.ref2querypos(read, x), range(vcf_record.start, vcf_record.stop)))
                    except IndexError:
                        read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
                    else:
                        if (mut_rng[0][1] != c.Ops.MATCH.value or
                            mut_rng[-1][1] != c.Ops.MATCH.value or
                            any(x[1] != c.Ops.DEL.value for x in mut_rng[1:-2])):
                            read_flag |= c.ValidatorFlags.BAD_OP.value
                elif vcf_record.rlen == 1:  # INS
                    try:
                        mut_rng = list(map(lambda x: r2s.ref2querypos(read, x), range(vcf_record.start, (vcf_record.start + len(alt)))))
                    except IndexError:
                        read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
                    else:
                        if (mut_rng[0][1] != c.Ops.MATCH.value or
                            mut_rng[-1][1] != c.Ops.MATCH.value or
                            any(x[1] != c.Ops.INS.value for x in mut_rng[1:-2])):
                            read_flag |= c.ValidatorFlags.BAD_OP.value
                        if read.query_sequence[mut_pos:len(alt)] != alt:  # type: ignore
                            read_flag |= c.ValidatorFlags.NOT_ALT.value
                else:  # COMPLEX
                    max_rng = range(vcf_record.start, vcf_record.stop) if (vcf_record.start + vcf_record.rlen) > (vcf_record.start + len(alt)) else range(vcf_record.start, (vcf_record.start + len(alt)))
                    try:
                        mut_rng = list(map(lambda x: r2s.ref2querypos(read, x), max_rng))
                    except IndexError:
                        read_flag |= c.ValidatorFlags.NOT_ALIGNED.value
                    else:
                        if (mut_rng[0][1] != c.Ops.MATCH.value or
                            mut_rng[-1][1] != c.Ops.MATCH.value):
                            read_flag |= c.ValidatorFlags.BAD_OP.value
                        if read.query_sequence[mut_pos:len(alt)] != alt:  # type: ignore
                            read_flag |= c.ValidatorFlags.NOT_ALT.value

                # n.b. nothing done if complex read
        if read_flag == c.ValidatorFlags.CLEAR.value:
            # "next", through an unfortunate quirk of history, means "mate", so this is reliable (pulls RNEXT)
            mate_end = r2s.ref_end_via_cigar(mate_cig, read.next_reference_start)  # type:ignore
            if not (read.flag & 0x40):
                # this looks like it should be checked for indexing snags
                pair_start = read.reference_start
                pair_end = read.reference_end
                if read.flag & 0x10:
                    if pair_start <= mate_end:
                        pair_start = mate_end + 1
                else:
                    if pair_end >= read.next_reference_start:  # type:ignore
                        pair_end = read.next_reference_start - 1
                if not (pair_start <= vcf_record.start <= pair_end):  # type:ignore
                    read_flag |= c.ValidatorFlags.OVERLAP.value
    return read_flag


def test_variant(
    vcf_rec: pysam.VariantRecord,
    mutant_bams: dict[str, pysam.AlignmentFile],
    alt: str,
    al_thresh: float,
    max_span: int,
    position_fraction_thresh: float,
    read_validator: c.FlagReturn,
) -> c.Filters:

    hp_filt = c.HPFilter()
    al_filt = c.ALFilter()

    mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in mutant_bams}
    mut_reads_log: dict[str, list[tuple]] = {key: [] for key in mutant_bams}
    mut_read_pos_f: list[int] = []
    mut_read_pos_r: list[int] = []
    mut_read_fracs_f: list[float] = []
    mut_read_fracs_r: list[float] = []
    aln_scores: list[float] = []

    for mut_sample, bam in mutant_bams.items():
        read_iter, test_iter = tee(bam.fetch(vcf_rec.chrom, vcf_rec.start, (vcf_rec.start + 1)))
        try:
            next(test_iter)
        except StopIteration:
            continue
        sample_readpair_ends = []
        read = None
        for read in read_iter: # type: ignore
            read_flag = c.ValidatorFlags.CLEAR.value
            read_flag = read_validator(vcf_record=vcf_rec, read=read, alt=alt)

            if read_flag == c.ValidatorFlags.CLEAR.value:
                mut_reads[mut_sample].append(read)
                sample_readpair_ends.append([read.reference_start, read.reference_end, read.next_reference_start, r2s.ref_end_via_cigar(read.get_tag('MC'), read.next_reference_start)])  # type: ignore
            mut_reads_log[mut_sample].append((read.query_name, read_flag))
            del(read)
        if len(mut_reads[mut_sample]) > 1:
            sample_readpair_ends_sorted: list[list[int]] = sorted(list(map(sorted, sample_readpair_ends)))
            curr_ends = [sample_readpair_ends_sorted[0]]
            drop_idx = []
            for i in range(1, len(sample_readpair_ends_sorted)):
                max_spans = map(lambda sublist: max([abs(x - y) for x, y in zip(sublist, sample_readpair_ends_sorted[i])]), curr_ends)
                if all([x <= max_span for x in max_spans]):
                    curr_ends.append(sample_readpair_ends_sorted[i])
                    drop_idx.append(i)
                else:
                    curr_ends = [sample_readpair_ends_sorted[i]]
            mut_reads[mut_sample] = [j for i, j in enumerate(mut_reads[mut_sample]) if i not in drop_idx]
    if all([len(x) == 0 for x in mut_reads.values()]):
        al_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
        hp_filt.code = c.FiltCodes.INSUFFICIENT_READS.value
    else:
        for read_list in mut_reads.values():
            for read in read_list:
                mut_pos, _ = r2s.ref2querypos(read, vcf_rec.start)
                if read.flag & 0x10:
                    read_idx_wrt_aln = read.query_alignment_end - mut_pos  # 1-based position where start, idx 1, is alignment end
                    mut_read_fracs_r.append(read_idx_wrt_aln / read.query_alignment_length)
                    mut_read_pos_r.append(read_idx_wrt_aln)
                else:
                    read_idx_wrt_aln  = mut_pos - read.query_alignment_start + 1
                    mut_read_fracs_f.append(read_idx_wrt_aln / read.query_alignment_length)
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
            if (((sum([x <= position_fraction_thresh for x in mut_read_fracs_f]) / len(mut_read_pos_f)) < 0.9) and
                  mad_f > 0 and
                  sd_f > 4):
                hp_filt.code = c.FiltCodes.SIXTYAI.value  # 60A(i)
            else:
                hp_filt.code = c.FiltCodes.SIXTYAI.value
                hp_filt.set()
        elif len(mut_read_pos_r) > 1 and not len(mut_read_pos_f) > 1:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            if (((sum([x <= position_fraction_thresh for x in mut_read_fracs_r]) / len(mut_read_pos_r)) < 0.9) and
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
            frac_lt_thresh = sum([x <= position_fraction_thresh for x in mut_read_fracs_f + mut_read_fracs_r]) / (len(mut_read_pos_f) + len(mut_read_pos_r))
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
    bams: dict[str, pysam.AlignmentFile],
    vcf_rec: pysam.VariantRecord,
    variant_tester: c.FiltReturn,
) -> dict[str, c.Filters]:

    if vcf_rec.alts is None:
        raise c.NoAlts
    samples_w_mutants = [name for name in vcf_rec.samples if vcf_rec.samples[name]["GT"] != (0, 0)]
    if len(samples_w_mutants) == 0:
        raise c.NoMutants

    bams_w_mutants = {k: v for k, v in bams.items() if k in samples_w_mutants}
    filt_d = {}
    for alt in vcf_rec.alts:
        filt_d[alt] = variant_tester(vcf_rec, bams_w_mutants, alt)
    return filt_d


def main_cli() -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s', datefmt='%I:%M:%S')

    parser = argparse.ArgumentParser(prog="hairpin")
    parser._optionals.title = 'info'
    parser.add_argument('-v', '--version', help='print version', action='version', version=c.VERSION)
    req = parser.add_argument_group('basic')
    req.add_argument('-i', '--vcf-in', help="path to input vcf")
    req.add_argument('-o', '--vcf-out', help="path to vcf out")
    req.add_argument('-b', '--bams', help="list of paths to name-sorted bams for samples in input vcf, whitespace separated", nargs='+')
    opt = parser.add_argument_group('extended')
    opt.add_argument('-ji', '--input-json', help='path to JSON of input parameters; overridden by arguments provided on command line', type=str)
    opt.add_argument('-jo', '--output-json', help='log input arguments to JSON', type=str)
    opt.add_argument('-m', '--name-mapping', help='map VCF sample names to BAM sample names; useful if they differ', metavar='VCF:BAM', nargs='+')
    opt.add_argument('-al', '--al-filter-threshold', help='threshhold for median of read alignment scores over read length, below which a variant is flagged as ALF - default: 0.93', type=float)
    opt.add_argument('-mc', '--min-clip-quality', help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35', type=int)
    opt.add_argument('-mq', '--min-mapping-quality', help='discard reads with mapping quality below this value - default: 11', type=int)
    opt.add_argument('-mb', '--min-base-quality', help='discard reads with base quality at variant position below this value - default: 25', type=int )
    opt.add_argument('-ms', '--max-read-span', help='maximum +- position to use when detecting PCR duplicates - default: 6', type=int)
    opt.add_argument('-pf', '--position-fraction', help='>90%% of variant reads variant must occur within [fraction] of start/end to allow HPF flag - default: 0.15', type=float)

    args = parser.parse_args()

    json_config: dict | None = None
    if args.input_json:
        logging.info('args JSON provided, arguments will be loaded from JSON if not present on command line')
        try:
            with open(args.input_json, 'r') as f:
                json_config = json.load(f)
        except Exception as e:
            h.cleanup(msg='failed to open input JSON, reporting: {}'.format(e))
        else:
            if not h.verify_json(json_config): # type:ignore
                h.cleanup(msg='JSON keys are not subset of available arguments (excluding --input-json and --output_json)')

    # set arg defaults
    for k in vars(args).keys():
        if not vars(args)[k]:
            if json_config and k in json_config.keys():
                setattr(args, k, json_config[k])
            elif k in c.DEFAULTS.keys():
                setattr(args, k, c.DEFAULTS[k])

    # test args are sensible, exit if not
    h.test_options(args)

    if args.output_json:
        try:
            with open(args.output_json, "w") as output_json:
                json.dump({k: vars(args)[k] for k in (vars(args).keys() - {'input_json', 'output_json'})}, output_json, indent="")
        except Exception as e:
            h.cleanup(msg='failed to write output JSON, reporting: {}'.format(e))

    primed_validate_read = partial(validate_read,
                                   min_mapqual=args.min_mapping_quality,
                                   min_clipqual=args.min_clip_quality,
                                   min_basequal=args.min_base_quality)

    primed_variant_tester = partial(test_variant, al_thresh=args.al_filter_threshold, max_span=args.max_read_span, position_fraction_thresh=args.position_fraction, read_validator=primed_validate_read)

    try:
        vcf_in_handle = pysam.VariantFile(args.vcf_in)
    except Exception as e:
        h.cleanup(msg='failed to open VCF input, reporting: {}'.format(e))

    sample_names = list(vcf_in_handle.header.samples)  # type:ignore
    if len(set(sample_names)) != len(sample_names):
        h.cleanup(msg='duplicate sample names in VCF')
    sample_names: set[str] = set(sample_names)
    vcf_sample_to_bam_file: dict[str, pysam.AlignmentFile] = {}
    for path in args.bams:
        try:
            bam = pysam.AlignmentFile(path, 'rb')
        except Exception as e:
            h.cleanup(msg='failed to read BAM at {}, reporting: {}'.format(path, e))
        # grab the sample name from first SM field
        # in header field RG
        # this may cause problems?
        # check with Peter
        bam_sample_name = bam.header.to_dict()['RG'][0]['SM']  # type:ignore
        vcf_sample_to_bam_file[bam_sample_name] = bam  # type:ignore
    if args.name_mapping:
        if len(args.name_mapping) > len(args.bams):
            h.cleanup(msg="more name mappings provided than BAMs")
        vcf_map_names = []
        bam_map_names = []
        for pair in args.name_mapping:
            kv_split = pair.split(':')  # VCF:BAM
            if len(kv_split) != 2:
                h.cleanup(msg='name mapping misformatted, more than two elements in map string {}'.format(pair))
            vcf_map_names.append(kv_split[0])
            bam_map_names.append(kv_split[1])
        if h.has_duplicates(vcf_map_names):
            h.cleanup(msg='duplicate VCF sample names provided to name mapping flag')
        if not set(vcf_map_names) <= sample_names:
            h.cleanup(msg="VCF sample names provided to name mapping flag are not equal to, or a subset of, VCF sample names as retrieved from VCF")
        if h.has_duplicates(bam_map_names):
            h.cleanup(msg='duplicate BAM sample names provided to name mapping flag')
        if h.lists_not_equal(bam_map_names, vcf_sample_to_bam_file.keys()):  # type:ignore
            h.cleanup(msg='BAM sample names provided to name mapping flag do not match BAM SM tags')
        vcf_sample_to_bam_file = {vcf_map_names[bam_map_names.index(k)]: v for k, v in vcf_sample_to_bam_file.items()}
    else:
        if not vcf_sample_to_bam_file.keys() <= sample_names:
            h.cleanup(msg='BAM SM tags do not match VCF sample names: {}'.format(vcf_sample_to_bam_file.keys() - sample_names))
    if sample_names != vcf_sample_to_bam_file.keys():
        logging.info("BAMs not provided for all VCF samples; {} will be ignored".format(sample_names - vcf_sample_to_bam_file.keys()))

    # init output
    out_head = vcf_in_handle.header  # type:ignore
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}, using samples {}\">".format(args.al_filter_threshold, ', '.join(vcf_sample_to_bam_file.keys())))
    out_head.add_line("##FILTER=<ID=HPF,Description=\"Variant arises from hairpin artefact, using samples {}\">".format(', '.join(vcf_sample_to_bam_file.keys())))
    out_head.add_line("##INFO=<ID=HPF,Number=1,Type=String,Description=\"alt|code for each alt indicating hairpin filter decision code\">")
    out_head.add_line("##INFO=<ID=ALF,Number=1,Type=String,Description=\"alt|code|score for each alt indicating AL filter conditions\">")

    try:
        vcf_out_handle = pysam.VariantFile(args.vcf_out, 'w', header=out_head)
    except Exception as e:
        h.cleanup(msg='failed to open VCF output, reporting: {}'.format(e))

    for record in vcf_in_handle.fetch():  # type:ignore
        try:
            filter_d: dict[str, c.Filters] = test_record_per_alt(
                bams=vcf_sample_to_bam_file,
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
                    record.info.update({filter.name: '|'.join([alt] + [str(f) if not type(f) == float else str(round(f, 3)) for f in filter][2:])})

            try:
                vcf_out_handle.write(record)  # type:ignore
            except Exception as e:
                h.cleanup(msg='failed to write to vcf, reporting: {}'.format(e))
    h.cleanup(c.EXIT_SUCCESS)
