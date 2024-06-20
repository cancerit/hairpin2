import pysam
from enum import Enum
from hairpin import ref2seq
from statistics import mean, median, stdev
import argparse
import logging
from itertools import tee
import sys

Ops = Enum('Ops',
           ['match', 'ins', 'delete', 'skip', 'soft', 'hard', 'pad', 'equal', 'diff', 'back'],
           start = 0)


def main(
    bam_paths: list,
    vcf_in_path: str,
    vcf_out_path: str,
    clip_qual_cutoff: int,
    min_mapqual: int,
    min_basequal: int,
    max_span: int,
    al_thresh: float,
    cent90_thresh:float
) -> None:
    
    vcf_obj = pysam.VariantFile(vcf_in_path)
    # init output
    out_head = vcf_obj.header
    out_head.add_line("##FILTER=<ID=ALF,Description=\"Median alignment score of reads reporting variant less than {}\">".format(al_thresh))
    out_head.add_line("##FILTER=<ID=HPF,Description=\"Evidence that variant arises from hairpin artefact\">")
    vcf_out = pysam.VariantFile(vcf_out_path, 'w', header=out_head)
    
    sample_names: list[str] = list(vcf_obj.header.samples)

    # add try excepts
    
    # func for local scope
    bam_reader_dict: dict[str, pysam.AlignmentFile] = dict.fromkeys(sample_names)
    def init_bam_dict() -> None:
        for path in bam_paths:
            bam = pysam.AlignmentFile(path, 'rb')
            # grab the sample name from first SM field
            # in header field RG
            # this may cause problems?
            # check with Peter
            bam_sample = bam.header.to_dict()['RG'][1]['SM']
            if bam_sample not in sample_names:
                logging.error('bam doesnt match')
                exit(1) # error
            else:
                bam_reader_dict[bam_sample] = bam
        if any([x is None for x in bam_reader_dict.values()]):
            logging.error('not enough bams')
            exit(1)
    init_bam_dict()  # execute w/o polluting namespace 
    
    mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in sample_names}
    # since intervals are unnecessary
    # - they were an artifact of the shearwater mp -
    # just iterate through all records in the vcf
    for vcf_rec in vcf_obj.fetch():

        if vcf_rec.alts is None:
            logging.error('vcf rec has no alts')
            continue  # ? ask Peter

        alt_test: bool = len(vcf_rec.alts[0]) == 1
        if vcf_rec.rlen == 1:
            mut_type = "sub" if alt_test else "ins"
        elif alt_test:
            mut_type = "del"
        else:
            mut_type = "complex"
            
        # check with Peter
        samples_w_mutants = [name for name in sample_names if vcf_rec.samples[name]["GT"] == (0, 1)]

        if len(samples_w_mutants) == 0:
            logging.error('no mutants')
            sys.exit(1)
        for mut_sample_name in samples_w_mutants:
            ### get_mutant_reads
            read_iter, test = tee(bam_reader_dict[mut_sample_name].fetch(vcf_rec.chrom, (vcf_rec.pos - 2), vcf_rec.pos))
            try:
                next(test)
            except StopAsyncIteration:
                logging.error('empty iterator')
                sys.exit(1)
            for read in read_iter: # type: ignore
                if any(x is None for x in [read.query_sequence, read.query_qualities, read.cigarstring, read.reference_start, read.reference_end]):
                    # breakpoint()
                    continue  # ?
            
                if not (read.flag & 0x2) or read.flag & 0xE00 or read.mapping_quality < min_mapqual:
                    # breakpoint()
                    continue
            
                mut_pos = ref2seq.ref2querypos(read, (vcf_rec.pos - 1)) # VCF 1-INDEXED, BAM 0-INDEXED
                mut_op = ref2seq.pos2op(mut_pos, read) if mut_pos != -1 else None

                # Check whether read reports variant or not - skip the rest of the loop if does not report variant
                # First, check for sub
                # does read.query_x work? or should it be read.query_alignment_x?
                if (mut_type == "sub" and
                    (not (mut_op == Ops.match.value or mut_op == Ops.diff.value) or
                     read.query_sequence[mut_pos] != vcf_rec.alts[0] or
                     read.query_qualities[mut_pos] < min_basequal)):
                    # breakpoint()
                    continue
                # Second, check whether length of read can accommodate size of indel
                # what if other alt is longer?
                if (mut_pos + vcf_rec.rlen > read.query_length or
                    mut_pos + len(vcf_rec.alts[0]) > read.query_length):
                    # breakpoint()
                    continue

                if mut_type == "del":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range((vcf_rec.pos - 1), vcf_rec.pos))  # check with peter re in/exclusivity of range
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] != Ops.match.value or
                        mut_rng_ops[-1] != Ops.match.value or
                        any(x != Ops.delete.value for x in mut_rng_ops)):
                        # breakpoint()
                        continue
                elif mut_type == "ins":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range(mut_pos, mut_pos + len(vcf_rec.alts[0]) + 1))
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] != Ops.match.value or
                        mut_rng_ops[-1] != Ops.match.value or
                        any(x != Ops.ins.value for x in mut_rng_ops) or
                        read.query_sequence[mut_pos:len(vcf_rec.alts[0])] != vcf_rec.alts[0]):
                        # breakpoint()
                        continue
            
                if ('S' in read.cigarstring and  # type: ignore
                    mean(read.query_alignment_qualities) < clip_qual_cutoff):  # type: ignore
                    # breakpoint()
                    continue

                if read.flag & 0x40:  # read first in pair
                    mut_reads[mut_sample_name].append(read)
                else:  # read second in pair
                    read_start = read.reference_start
                    read_end = read.reference_end
                    if read.flag & 0x10:
                        mate = bam_reader_dict[mut_sample_name].mate(read)
                        if mate.reference_end is None:
                            # breakpoint()
                            continue  # ?
                        if read.reference_start <= mate.reference_end:
                            read_start = mate.reference_end + 1
                    else:
                        if read.reference_end >= read.next_reference_start:  # check with Peter, does this map on to his code accurately
                            read_end = read.next_reference_start - 1
                    # Peter's implemenation comments that this conditional below
                    # checks if mutant overlaps with first read in pair
                    # I don't see it, perhaps it refers to the code above
                    if read_start <= (vcf_rec.pos - 1) <= read_end:
                        mut_reads[mut_sample_name].append(read)
                    else:
                        # breakpoint() 
                        continue
        del(read_iter, test, read, mut_sample_name)  # type: ignore
        # breakpoint()
        if all([x is None for x in mut_reads.values()]):
            logging.error('empty mut_reads')
            sys.exit(1)
        
        ### remove_dups_with_wobble
        for mut_sample_name, reads in mut_reads.items():
            if len(reads) == 0:
                continue
            # want the start of the record, the end
            ### start_mate_end_pairs()
            # incidentally, I suppose hairpin only works for paired data?
            sorted_ends = []
            for read in reads:
                mate = bam_reader_dict[mut_sample_name].mate(read)
            
                if any(x is None for x in [read.reference_start,
                                        read.reference_end,
                                        read.cigartuples,
                                        mate.reference_start,
                                        mate.reference_end,
                                        mate.cigartuples]):
                    continue  # ?
            
                # this gets pos wrt to alignment against reference
                start: int = read.reference_start
                end: int = read.reference_end
                mate_start: int  = mate.reference_start
                mate_end: int  = mate.reference_end
            
                # Peter
                # behaviour on cig = None?     
                # Julia XAM gets cigar info differently, checking CG:B,I tag
                # does this matter?
                if read.cigartuples[0][0] == Ops.soft.value:
                    start -= read.cigartuples[0][1]
                if read.cigartuples[-1][0] == Ops.soft.value:
                    end += read.cigartuples[-1][1]
                if mate.cigartuples[0][0] == Ops.soft.value:
                    mate_start -= mate.cigartuples[0][1]
                if mate.cigartuples[-1][0] == Ops.soft.value:
                    mate_end += mate.cigartuples[-1][1]
            
                # appears mate posns simply aren't assigned if none
                sorted_ends.append(sorted([start, end, mate_start, mate_end]))
            ### end
            sorted_ends: list[list[int]] = sorted(sorted_ends)  # sort sublists on first element
            min_ends: list[list[int]] = [sorted_ends.pop(0)]
            # I dont' fully understand this segment but I think it recapitulates the Julia (ask Peter)
            i = 1
            while len(sorted_ends) != 0:
                loop_ends: list[int] = sorted_ends.pop(0)
                max_spans = map(lambda sublist: max([abs(x - y) for x, y in zip(sublist, loop_ends)]), min_ends)
             
                if all([x <= max_span for x in max_spans]):
                    min_ends.append(loop_ends)
                    reads.pop(i)
                else:
                    min_ends = [loop_ends]
                i += 1
            mut_reads[mut_sample_name] = reads
        del(read, reads, mut_sample_name)
        # breakpoint()
    
        ### check hairpin filter
        mut_read_pos_f: list[int] = []
        mut_read_pos_r: list[int] = []
        mut_read_fracs_f: list[float] = []
        mut_read_fracs_r: list[float] = []
        aln_scores: list[float] = []
        # for my test vcf, by the time we get here there is only on read, for one sample, in mut_reads
        # one sample is correct, since 0/1 only appears on in l0045
        # but Peter's implementation flags HPF, so perhaps there should be more samples?
        # need to follow that bam through this process, see what happens to each read
        for _, reads in mut_reads.items():
            for read in reads:
                # breakpoint()
                if any([x is None for x in [read.reference_start, read.reference_end]]):
                    # breakpoint()
                    continue  # ?
                
                mut_pos = ref2seq.ref2querypos(read, (vcf_rec.pos - 1))
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
                    aln_scores.append(read.get_tag('AS') / read.query_length)  # or should this be .query_alignment_length? (Peter)
                except KeyError:
                    # breakpoint()
                    continue  # ?
        al_filt = median(aln_scores) <= al_thresh

        # breakpoint()
        fbool = len(mut_read_pos_f) > 1
        rbool = len(mut_read_pos_r) > 1
        if fbool:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
        if rbool:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
        # hairpin conditions from Ellis et al.
        hp_filt = True
        # these branches all lead to the same result!
        if fbool and rbool:
            frac_lt_thresh = sum([x <= cent90_thresh for x in mut_read_fracs_f + mut_read_fracs_r]) / (len(mut_read_pos_f) + len(mut_read_pos_r))
            if (frac_lt_thresh < 0.9 or
                (mad_f > 2 and mad_r > 2 and sd_f > 2 and sd_r > 2) or
                (mad_f > 1 and sd_f > 10) or
                (mad_r > 1 and sd_r > 10)):
                # breakpoint()
                hp_filt = False
        elif fbool:
            if (((sum([x <= cent90_thresh for x in mut_read_pos_f]) / len(mut_read_pos_f)) < 0.9) and
                mad_f > 0 and
                sd_f > 4):
                # breakpoint()
                hp_filt = False
        elif rbool:
            if (((sum([x <= cent90_thresh for x in mut_read_pos_r]) / len(mut_read_pos_r)) < 0.9) and
                mad_r > 0 and
                sd_r > 4):
                # breakpoint()
                hp_filt = False
        else:
            # breakpoint()
            hp_filt = False
        ### end 

        ### update vcf record
        if al_filt:
            vcf_rec.filter.add("ALF")
        if hp_filt:
            vcf_rec.filter.add("HPF")
            breakpoint()
        
        # try except
        # breakpoint()
        vcf_out.write(vcf_rec)


if __name__ == '__main__':
    
    logging.basicConfig(level=logging.INFO)
    
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
    
    args = parser.parse_args()
    if any([x is None for _, x in vars(args).items()]):
        logging.info('option(s) not provided, using defaults')
    
    main(
        bam_paths=args.bams,
        vcf_in_path=args.vcf_in,
        vcf_out_path=args.vcf_out,
        clip_qual_cutoff=args.clip_quality_cutoff if args.clip_quality_cutoff else 35,
        min_mapqual=args.min_mapping_quality if args.min_mapping_quality else 11,
        min_basequal=args.min_base_quality if args.min_base_quality else 25,
        max_span=args.max_read_span if args.max_read_span else 6,
        al_thresh=args.AL_filter_threshold if args.AL_filter_threshold else 0.93,
        cent90_thresh=args.cent90_threshold if args.cent90_threshold else 0.15
    )
    
    
