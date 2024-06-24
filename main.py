import pysam
from enum import Enum
from hairpin import ref2seq, utilities
from statistics import mean, median, stdev
import argparse
import logging
from itertools import tee
from functools import partial
from sys import exit as sysexit
from dataclasses import dataclass
from typing import Callable, Optional

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

Ops = Enum('Ops',
           ['match', 'ins', 'delete', 'skip', 'soft', 'hard', 'pad', 'equal', 'diff', 'back'],
           start = 0)


# dataclass for tracking info?
# define exit codes and do all logging outside of main?


def cleanup(code: int, msg: Optional[str] = None) -> None:
    if code != EXIT_SUCCESS and msg:
        logging.error(msg)
    for obj_name in ['vcf_in_handle', 'vcf_out_handle']:
        if obj_name in locals():
            locals()[obj_name].close()  # lol
    if 'bam_reader_d' in locals():
        locals()['bam_reader_d'].close()
    if 'log_file' in locals() and locals()['log_file']:
        locals()['log_file'].close()
    if code == EXIT_SUCCESS:
        logging.info('hairpin complete')
    sysexit(code)


def main(
    bams: dict[str, pysam.AlignmentFile],
    vcf_in: pysam.VariantFile,
    vcf_out: pysam.VariantFile,
    clip_qual_cutoff: int,
    min_mapqual: int,
    min_basequal: int,
    max_span: int,
    al_thresh: float,
    cent90_thresh: float,
    log_func: Callable,
) -> int:
    
    for vcf_rec in vcf_in.fetch():

        if vcf_rec.alts is None:
            logging.error('Mutation {}:{} ¦ no alts in VCF')
            # Peter what is the correct approach to this state?
            continue  # ? 

        # Peter is it possible for other alts to be longer?
        # the vcf format doesn't specify as far as I can tell that
        # alts should be the same length (and elsewhere)
        if vcf_rec.rlen == 1:
            mut_type = "sub" if len(vcf_rec.alts[0]) == 1 else "ins"
        elif len(vcf_rec.alts[0]) == 1:
            mut_type = "del"
        else:
            mut_type = "complex"
            

        # Peter is this sufficiently specific? You search for 0/1 as a string
        samples_w_mutants = [name for name in sample_names if vcf_rec.samples[name]["GT"] == (0, 1)]

        if len(samples_w_mutants) == 0:
            logging.error('Mutation {}:{} --- no samples exhibiting mutation')
            # Peter what is the correct approach to this situtation?
            return EXIT_FAILURE

        mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in samples_w_mutants}

        mut_read_pos_f: list[int] = []
        mut_read_pos_r: list[int] = []
        mut_read_fracs_f: list[float] = []
        mut_read_fracs_r: list[float] = []
        aln_scores: list[float] = []
        
        for mut_sample_name in samples_w_mutants:
            read_iter, test_iter = tee(bams[mut_sample_name].fetch(vcf_rec.chrom, vcf_rec.start, (vcf_rec.start + 1)))
            try:
                next(test_iter)
            except StopIteration:
                logging.error('Mutation {}:{} --- no reads mapped to region despite 0/1 for sample {}'.format(vcf_rec.chrom, vcf_rec.pos, mut_sample_name))
                continue
            # Peter the above will skip the sample if no reads are returned by the fetch
            # is that acceptable? Would we not expect reads in the bam if the sample
            # exhibits 0/1
            sorted_ends = []
            read = None
            for read in read_iter: # type: ignore

                if not (read.flag & 0x2) or read.flag & 0xE00 or read.mapping_quality < min_mapqual:
                    continue

                if any(x is None for x in [read.reference_start,
                                        read.reference_end,
                                        read.query_sequence,
                                        read.query_qualities,
                                        read.cigarstring,
                                        read.cigartuples]):
                                        # Peter what is the correct response to this state?
                                        continue

                mut_pos = ref2seq.ref2querypos(read, vcf_rec.start) # VCF 1-INDEXED, BAM 0-INDEXED
                # Peter this occurs when the cigar string for a read
                # indicates a deletion over the reference position in the read
                # despite the vcf calling a single-base substitution
                # What is the correct approach to this state?
                if mut_pos is None:
                    continue
                mut_op = ref2seq.pos2op(mut_pos, read) if mut_pos != -1 else None


                # Check whether read reports variant or not - skip the rest of the loop if does not report variant
                # First, check for sub
                if mut_type == "sub":
                    if (mut_op not in [Ops.match.value, Ops.diff.value] or
                        read.query_sequence[mut_pos] != vcf_rec.alts[0] or  # what about other alts?
                        read.query_qualities[mut_pos] < min_basequal):
                            # breakpoint()
                            continue
                # Second, check whether length of read can accommodate size of indel
                # Peter again could other alt is longer?
                elif (mut_pos + vcf_rec.rlen > read.query_length or
                      mut_pos + len(vcf_rec.alts[0]) > read.query_length):
                    # breakpoint()
                    continue

                if mut_type == "del":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range(vcf_rec.start, vcf_rec.stop))
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] not in [Ops.match.value, Ops.diff.value] or  # Peter should diff op be considered here?
                        mut_rng_ops[-1] not in [Ops.match.value, Ops.diff.value] or
                        any(x != Ops.delete.value for x in mut_rng_ops[1:-2])):
                        # breakpoint()
                        continue
                elif mut_type == "ins":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range(vcf_rec.start, (vcf_rec.start + len(vcf_rec.alts[0]))))  # other alt?
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] not in [Ops.match.value, Ops.diff.value] or
                        mut_rng_ops[-1] not in [Ops.match.value, Ops.diff.value] or
                        any(x != Ops.ins.value for x in mut_rng_ops[1:-2]) or
                        read.query_sequence[mut_pos:len(vcf_rec.alts[0])] != vcf_rec.alts[0]):  # further checks needed, should it be query_alignment_sequence?
                        # breakpoint()
                        continue

                # n.b. nothing done if complex read
                
                if ('S' in read.cigarstring and  # type: ignore
                    mean(read.query_alignment_qualities) < clip_qual_cutoff):  # type: ignore
                    # breakpoint()
                    continue

                mate = bams[mut_sample_name].mate(read)

                if any(x is None for x in [mate.reference_start,
                                            mate.reference_end,
                                            mate.cigartuples]):
                                        # Peter what is the correct approach to this state
                                        continue

                if read.flag & 0x40:  # read first in pair
                    mut_reads[mut_sample_name].append(read)
                else:  # read second in pair
                    read_start = read.reference_start
                    read_end = read.reference_end
                    if read.flag & 0x10:
                        if mate.reference_end is None:
                            # breakpoint()
                            continue  # ?
                        if read.reference_start <= mate.reference_end:
                            read_start = mate.reference_end + 1
                    else:
                        if read.reference_end >= read.next_reference_start:
                            read_end = read.next_reference_start - 1
                    if read_start <= vcf_rec.start <= read_end:
                        mut_reads[mut_sample_name].append(read)
                    else:
                        # breakpoint() 
                        continue
                # if we've got this far, then read has been added to mut_reads
                # get pos wrt to aligned region
                # Peter
                # behaviour on cig = None?     
                # Also, Julia XAM gets cigar info differently, checking CG:B,I tag
                # does this matter?
                soft_start = (read.reference_start - read.cigartuples[0][1]) if read.cigartuples[0][0] == Ops.soft.value else read.reference_start
                soft_end = (read.reference_end + read.cigartuples[-1][1]) if read.cigartuples[-1][0] == Ops.soft.value else read.reference_end
                soft_mate_start = (mate.reference_start - mate.cigartuples[0][1]) if mate.cigartuples[0][0] == Ops.soft.value else mate.reference_start
                soft_mate_end = (mate.reference_end + mate.cigartuples[-1][1]) if mate.cigartuples[-1][0] == Ops.soft.value else mate.reference_end

                sorted_ends.append(sorted([soft_start, soft_end, soft_mate_start, soft_mate_end]))
                
            # if we got nothing for that sample after cycling through all reads
            if len(mut_reads[mut_sample_name]) == 0:
                # Peter is this an error state? if the sample showed 0/1 would we expect viable reads supporting that mutation?
                continue

            # return to per sample
            sorted_ends: list[list[int]] = sorted(sorted_ends)  # sort sublists on first element
            # Peter what if sorted_ends contains only 1 item?
            min_ends: list[list[int]] = [sorted_ends.pop(0)]
            # Peter I dont' fully understand this segment but I think it recapitulates the Julia
            i = 1
            drop_idx = []
            while len(sorted_ends) != 0:
                loop_ends: list[int] = sorted_ends.pop(0)
                max_spans = map(lambda sublist: max([abs(x - y) for x, y in zip(sublist, loop_ends)]), min_ends)

                if all([x <= max_span for x in max_spans]):
                    min_ends.append(loop_ends)
                    drop_idx.append(i)
                else:
                    min_ends = [loop_ends]
                i += 1
            mut_reads[mut_sample_name] = [j for i, j in enumerate(mut_reads[mut_sample_name]) if i not in drop_idx]
            del(i, drop_idx)
            
            if read:
                del(read)

        # if we got nothing for all samples exhibiting 0/1 for this mutation
        # skip this mutation
        # Peter please confirm correct approach to this state
        if all([len(x) == 0 for x in mut_reads.values()]):
            continue

        for read_list in mut_reads.values():
            for read in read_list:
                mut_pos = ref2seq.ref2querypos(read, vcf_rec.start)
                if mut_pos == -1:
                    # breakpoint()
                    continue  # ?
                if read.flag & 0x10:
                    read_loc = read.reference_end - mut_pos + 1
                    mut_read_fracs_r.append(read_loc / (read.reference_start - read.reference_end + 1))
                    mut_read_pos_r.append(read_loc)
                else:
                    read_loc = (mut_pos - read.reference_start + 1)
                    mut_read_fracs_f.append(read_loc / (read.reference_end - read.reference_start + 1))
                    mut_read_pos_f.append(read_loc)
                # breakpoint()
                try:
                    aln_scores.append(read.get_tag('AS') / read.query_length)
                except KeyError:
                    # Peter what is the correct approach to this state?
                    continue  # ?
        if len(aln_scores) == 0:
            breakpoint()
        al_filt = (avg_AS := median(aln_scores) <= al_thresh)

        fbool = len(mut_read_pos_f) > 1
        rbool = len(mut_read_pos_r) > 1
        # hairpin conditions from Ellis et al.
        if fbool and not rbool:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
            if (((sum([x <= cent90_thresh for x in mut_read_pos_f]) / len(mut_read_pos_f)) < 0.9) and
                  mad_f > 0 and
                  sd_f > 4):
                log_func(msg='Mutation {}:{} --- passed HPF per Ellis et al. 60B(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=2)
                hp_filt = False
            else:
                log_func(msg='Mutation {}:{} --- failed HPF per Ellis et al. 60B(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=1)
                hp_filt = True
        elif rbool and not fbool:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            if (((sum([x <= cent90_thresh for x in mut_read_pos_r]) / len(mut_read_pos_r)) < 0.9) and
                  mad_r > 0 and
                  sd_r > 4):
                log_func(msg='Mutation {}:{} --- passed HPF per Ellis et al. 60A(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=2)
                hp_filt = False
            else:
                log_func(msg='Mutation {}:{} --- failed HPF per Ellis et al. 60(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=1)
                hp_filt = True
        elif fbool and rbool:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
            frac_lt_thresh = sum([x <= cent90_thresh for x in mut_read_fracs_f + mut_read_fracs_r]) / (len(mut_read_pos_f) + len(mut_read_pos_r))
            if (frac_lt_thresh < 0.9 or
               (mad_f > 2 and mad_r > 2 and sd_f > 2 and sd_r > 2) or
               (mad_f > 1 and sd_f > 10) or
               (mad_r > 1 and sd_r > 10)):
                log_func(msg='Mutation {}:{} --- passed HPF per Ellis et al. 60A(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=2)
                hp_filt = False
            else:
                log_func(msg='Mutation {}:{} --- failed HPF per Ellis et al. 60A(i)'.format(vcf_rec.chrom, vcf_rec.pos), decision_lvl=1)
                hp_filt = True
        else:
            log_func(msg='Mutation {}:{} --- passed HPF, insufficient reads to support HPF flag', decision_lvl=2)  # Peter does this comment accurately assses the situation
            hp_filt = False

        ### update vcf record
        if al_filt:
            log_func(msg='Mutation {}:{} --- failed ALF, median AS of {}'.format(vcf_rec.chrom, vcf_rec.pos, avg_AS), decision_lvl=1)
            vcf_rec.filter.add("ALF")
        if hp_filt:
            vcf_rec.filter.add("HPF")

        try:
            vcf_out.write(vcf_rec)
        except Exception as e:
            logging.error('Failed to write VCF, reporting: {}'.format(e))
            return EXIT_FAILURE

    return EXIT_SUCCESS


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s', datefmt='%I:%M:%S')

    parser = argparse.ArgumentParser(prog="hairpin")
    parser._optionals.title = 'info'
    parser.add_argument('-v', '--version', help='print version', action='version', version='hairpin 1.0.0')
    req = parser.add_argument_group('required')
    req.add_argument('-i', '--vcf-in', help="path to input vcf", required=True)
    req.add_argument('-o', '--vcf-out', help="path to vcf out", required=True)
    req.add_argument('-b', '--bams', help="list of paths to bams for samples in input vcf, whitespace separated", nargs='+', required=True)
    opt = parser.add_argument_group('options')
    opt.add_argument('-cq', '--clip-quality-cutoff', help='default: 35', type=int)
    opt.add_argument('-mq', '--min-mapping-quality', help='default: 11', type=int)
    opt.add_argument('-mb', '--min-base-quality', help='default: 25', type=int)
    opt.add_argument('-ms', '--max-read-span', help='default: 6', type=int)
    opt.add_argument('-al', '--AL-filter-threshold', help='default: 0.93', type=float)
    opt.add_argument('-c9', '--cent90-threshold', help='default: 0.15', type=float)
    log_sev_opt = opt.add_mutually_exclusive_group()
    log_sev_opt.add_argument('-l', dest='log_path', help='log reason for failing records to file', action=utilities.SetLogAndSeverity)
    log_sev_opt.add_argument('-ll', dest='log_path', help='as -l, and addtionally log reason for passing records')
    log_sev_opt.add_argument('-lll', dest='log_path', help='as -ll, and additionaly log reason for discarding reads associated with a record', action=utilities.SetLogAndSeverity)

    args = parser.parse_args()

    if not hasattr(args, 'severity'):
        args.severity = 0

    if any([x is None for _, x in vars(args).items()]):
        logging.info('option(s) not provided, using defaults')

    try:
        log_file = open(args.log_path) if args.log_path else None
    except Exception as e:
        cleanup(1, 'failed to open log file, reporting: {}'.format(e))
    primed_log_func = partial(utilities.log_decision, log_lvl=args.severity, log_file=log_file)

    try:
        vcf_in_handle = pysam.VariantFile(args.vcf_in)
    except Exception as e:
        cleanup(1, 'failed to open VCF input, reporting: {}'.format(e))

    # init output
    out_head = vcf_in_handle.header
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}\">".format(args.AL_filter_threshold if args.AL_filter_threshold else None))
    out_head.add_line("##FILTER=<ID=HPF,Description=\"Evidence that variant arises from hairpin artefact\">")

    try:
        vcf_out_handle = pysam.VariantFile(args.vcf_out, 'w', header=out_head)
    except Exception as e:
        cleanup(1, 'failed to open VCF output, reporting: {}'.format(e))

    sample_names: list[str] = list(vcf_out_handle.header.samples)

    # func for local scope
    bam_reader_d: dict[str, pysam.AlignmentFile] = dict.fromkeys(sample_names)
    for path in args.bams:
        try:
            bam = pysam.AlignmentFile(path, 'rb')
        except Exception as e:
            cleanup(1, 'failed to read BAM at {}, reporting: {}'.format(path, e))
        # grab the sample name from first SM field
        # in header field RG
        # this may cause problems?
        # check with Peter
        bam_sample = bam.header.to_dict()['RG'][1]['SM']
        if bam_sample not in sample_names:
            cleanup(1, 'name in header ({}) of BAM at {} does not match any samples in VCF'.format(bam_sample, path))
        else:
            bam_reader_d[bam_sample] = bam


    main_code = main(
        bams=bam_reader_d,
        vcf_in=vcf_in_handle,
        vcf_out=vcf_out_handle,
        clip_qual_cutoff=args.clip_quality_cutoff if args.clip_quality_cutoff else 35,
        min_mapqual=args.min_mapping_quality if args.min_mapping_quality else 11,
        min_basequal=args.min_base_quality if args.min_base_quality else 25,
        max_span=args.max_read_span if args.max_read_span else 6,
        al_thresh=args.AL_filter_threshold if args.AL_filter_threshold else 0.93,
        cent90_thresh=args.cent90_threshold if args.cent90_threshold else 0.15,
        log_func=primed_log_func
    )

    cleanup(main_code)
