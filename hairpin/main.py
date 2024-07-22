import pysam
from hairpin import ref2seq as r2s, constants as c
from statistics import mean, median, stdev
import argparse
import logging
import json
from itertools import tee
from functools import partial
import sys


def cleanup(code: int = c.EXIT_FAILURE, msg: None | str = None) -> None:
    if code != c.EXIT_SUCCESS and msg:
        logging.error(msg)
    for obj_name in ['vcf_in_handle', 'vcf_out_handle']:
        if obj_name in locals():
            locals()[obj_name].close()  # lol
    if 'bam_reader_d' in locals():
        locals()['bam_reader_d'].close()
    if 'log_file' in locals() and locals()['log_file']:
        locals()['log_file'].close()
    if code == c.EXIT_SUCCESS:
        logging.info('hairpin complete')
    sys.exit(code)


# CIGAR best retrieved from CG:B,I tag - implement in future
def validate_read(
    vcf_record: pysam.VariantRecord,
    read: pysam.AlignedSegment,
    min_mapqual: int,
    clip_qual_cutoff: int,
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
            mean(read.query_alignment_qualities) < clip_qual_cutoff):  # type: ignore
            read_flag |= c.ValidatorFlags.CLIPQUAL.value
        # First, check for sub
        try:
            mut_pos, mut_op = r2s.ref2querypos(read, vcf_record.start) # VCF 1-INDEXED, BAM 0-INDEXED - vcf_record.start = 0-indexed mutation position. testing with pos, 1-indexed, to see if match Peter
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
            # is it safe to assume this is always mate?
            mate_end = r2s.ref_end_via_cigar(mate_cig, read.next_reference_start)  # THIS ONLY WORKS ASSUMING MATE IS NEXT READ
            if not (read.flag & 0x40):
                # this looks like it should be checked for indexing snags
                pair_start = read.reference_start
                pair_end = read.reference_end
                if read.flag & 0x10:
                    if pair_start <= mate_end:
                        pair_start = mate_end + 1
                else:
                    if pair_end >= read.next_reference_start:
                        pair_end = read.next_reference_start - 1
                if not (pair_start <= vcf_record.start <= pair_end):
                    read_flag |= c.ValidatorFlags.OVERLAP.value
    return read_flag


def test_variant(
    vcf_rec: pysam.VariantRecord,
    mutant_bams: dict[str, pysam.AlignmentFile],
    alt: str,
    al_thresh: float,
    max_span: int,
    cent90_thresh: float,
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
                    aln_scores.append(read.get_tag('AS') / read.query_length)
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
            if (((sum([x <= cent90_thresh for x in mut_read_fracs_f]) / len(mut_read_pos_f)) < 0.9) and
                  mad_f > 0 and
                  sd_f > 4):
                hp_filt.code = c.FiltCodes.SIXTYAI.value  # 60A(i)
            else:
                hp_filt.code = c.FiltCodes.SIXTYAI.value
                hp_filt.set()
        elif len(mut_read_pos_r) > 1 and not len(mut_read_pos_f) > 1:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            if (((sum([x <= cent90_thresh for x in mut_read_fracs_r]) / len(mut_read_pos_r)) < 0.9) and
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
            frac_lt_thresh = sum([x <= cent90_thresh for x in mut_read_fracs_f + mut_read_fracs_r]) / (len(mut_read_pos_f) + len(mut_read_pos_r))
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
    req = parser.add_argument_group('required')
    req.add_argument('-i', '--vcf-in', help="path to input vcf", required=True)
    req.add_argument('-o', '--vcf-out', help="path to vcf out", required=True)
    req.add_argument('-b', '--bams', help="list of paths to bams for samples in input vcf, whitespace separated", nargs='+', required=True)
    opt = parser.add_argument_group('options')
    opt.add_argument('-cq', '--clip-quality-cutoff', help='default: 35', type=int, default=35)
    opt.add_argument('-mq', '--min-mapping-quality', help='default: 11', type=int, default=11)
    opt.add_argument('-mb', '--min-base-quality', help='default: 25', type=int, default=25)
    opt.add_argument('-ms', '--max-read-span', help='default: 6', type=int, default=6)
    opt.add_argument('-al', '--al-filter-threshold', help='default: 0.93', type=float, default=0.93)
    opt.add_argument('-c9', '--cent90-threshold', help='default: 0.15', type=float, default=0.15)
    opt.add_argument('-j', '--json-log', dest='json_path', help='log input parameters/arguments to JSON', type=str)

    args = parser.parse_args()

    primed_validate_read = partial(validate_read,
                                   min_mapqual=args.min_mapping_quality,
                                   clip_qual_cutoff=args.clip_quality_cutoff,
                                   min_basequal=args.min_base_quality)

    primed_variant_tester = partial(test_variant, al_thresh=args.al_filter_threshold, max_span=args.max_read_span, cent90_thresh=args.cent90_threshold, read_validator=primed_validate_read)

    try:
        vcf_in_handle = pysam.VariantFile(args.vcf_in)
    except Exception as e:
        cleanup(msg='failed to open VCF input, reporting: {}'.format(e))

    # init output
    out_head = vcf_in_handle.header
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}\">".format(args.al_filter_threshold))
    out_head.add_line("##FILTER=<ID=HPF,Description=\"Evidence that variant arises from hairpin artefact\">")
    out_head.add_line("##INFO=<ID=HPF,Number=1,Type=String,Description=\"alt|code for each alt indicating hairpin filter decision code\">")
    out_head.add_line("##INFO=<ID=ALF,Number=1,Type=String,Description=\"alt|code|score for each alt indicating AL filter conditions\">")

    try:
        vcf_out_handle = pysam.VariantFile(args.vcf_out, 'w', header=out_head)
    except Exception as e:
        cleanup(msg='failed to open VCF output, reporting: {}'.format(e))

    sample_names: list[str] = list(vcf_in_handle.header.samples)

    bam_reader_d: dict[str, None | pysam.AlignmentFile] = dict.fromkeys(sample_names)
    for path in args.bams:
        try:
            bam = pysam.AlignmentFile(path, 'rb')
        except Exception as e:
            cleanup(msg='failed to read BAM at {}, reporting: {}'.format(path, e))
        # grab the sample name from first SM field
        # in header field RG
        # this may cause problems?
        # check with Peter
        bam_sample = bam.header.to_dict()['RG'][1]['SM']
        if bam_sample not in sample_names:
            cleanup(msg='name in header ({}) of BAM at {} does not match any samples in VCF'.format(bam_sample, path))
        else:
            bam_reader_d[bam_sample] = bam

    for record in vcf_in_handle.fetch():
        try:
            filter_d: dict[str, c.Filters] = test_record_per_alt(
                bams=bam_reader_d,  # type: ignore
                vcf_rec=record,
                variant_tester=primed_variant_tester
            )
        except c.NoAlts:
            logging.warning('{0: <7}:{1: >12} ¦ no alts for this record'.format(record.chrom, record.pos))
        except c.NoMutants:
            logging.warning('{0: <7}:{1: >12} ¦ no samples contain reads exhibiting record alts'.format(record.chrom, record.pos))
        else:
            for alt, filter_bundle in filter_d.items():
                for filter in filter_bundle:
                    if filter.flag:
                        record.filter.add(filter.name)
                    record.info.update({filter.name: '|'.join([alt] + [str(f) if not type(f) == float else str(round(f, 3)) for f in filter][2:])})

            try:
                vcf_out_handle.write(record)
            except Exception as e:
                cleanup(msg='failed to write to vcf, reporting: {}'.format(e))
    if args.json_path:
        try:
            with open(args.json_path, "w") as jo:
                json.dump(vars(args), jo)
        except Exception as e:
            logging.warning('retaining output, but failed to write to parameters json, reporting {}'.format(e))
    cleanup(c.EXIT_SUCCESS)
