import pysam
from enum import Enum
from typing import Optional
from hairpin import ref2seq
from statistics import mean, median, stdev
import argparse

Ops = Enum('Ops',
           ['match', 'ins', 'delete', 'skip', 'soft', 'hard', 'pad', 'equal', 'diff', 'back'],
           start = 0)


# is streaming approach necessary?
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

    # try excepts
    bam_reader_dict: dict[str, Optional[pysam.AlignmentFile]] = dict.fromkeys(sample_names)
    for path in bam_paths:
        bam = pysam.AlignmentFile(path, 'rb')
        # grab the sample name from first SM field
        # in header field RG
        # this may cause problems?
        # check with Peter
        if bam_sample := bam.header.to_dict()['RG'][1]['SM'] not in sample_names:
            exit(1) # error
        else:
            bam_reader_dict[bam_sample] = bam  # type: ignore

    mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in sample_names}
    # since intervals are unnecessary
    # - they were an artifact of the shearwater mp -
    # just iterate through all records in the vcf
    for vcf_rec in vcf_obj.fetch():
        
        if vcf_rec.alts is None:
            continue  # ?
        alt_test: bool = len(vcf_rec.alts[0]) == 1
        if vcf_rec.rlen == 1:
            mut_type = "sub" if alt_test else "ins"
        elif alt_test:
            mut_type = "del"
        else:
            mut_type = "complex"
                
        # check with Peter
        samples_w_mutants = [name for name in sample_names if vcf_rec.samples[name]["GT"] == (0, 1)]
        
        for mut_sample_name in samples_w_mutants:
            ### get_mutant_reads
            for read in bam_reader_dict[mut_sample_name].fetch(vcf_rec.chrom, vcf_rec.start, vcf_rec.stop): # type: ignore
                
                if any(x is None for x in [read.query_sequence, read.query_qualities, read.cigarstring, read.reference_start, read.reference_end]):
                    continue  # ?
                
                if read.flag & 0xE02 or read.mapping_quality < min_mapqual:
                    continue
                
                mut_pos = ref2seq.ref2querypos(read, vcf_rec.pos)
                mut_op = ref2seq.pos2op(mut_pos, read) if mut_pos != -1 else None

                # Check whether read reports variant or not - skip the rest of the loop if does not report variant
                # First, check for sub
                # does read.query_x work? or should it be read.query_alignment_x?
                if (mut_type == "sub" and
                    (mut_op != Ops.match or
                     read.query_sequence[mut_pos] != vcf_rec.alts[0] or
                     read.query_qualities[mut_pos] < min_basequal)):
                    continue
                
                # Second, check whether length of read can accommodate size of indel
                # what if other alt is longer?
                if (mut_pos + vcf_rec.rlen > read.query_length or
                    mut_pos + len(vcf_rec.alts[0]) > read.query_length):
                    continue
                
                if mut_type == "del":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range(vcf_rec.start, vcf_rec.stop + 1))  # check with peter re in/exclusivity of range
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] != Ops.match or
                        mut_rng_ops[-1] != Ops.match or
                        any(x != Ops.delete for x in mut_rng_ops)):
                        continue
                elif mut_type == "ins":
                    mut_rng = map(lambda x: ref2seq.ref2querypos(read, x), range(mut_pos, mut_pos + len(vcf_rec.alts[0]) + 1))
                    mut_rng_ops = list(map(lambda x: ref2seq.pos2op(x, read), mut_rng))
                    if (mut_rng_ops[0] != Ops.match or
                        mut_rng_ops[-1] != Ops.match or
                        any(x != Ops.ins for x in mut_rng_ops) or
                        read.query_sequence[mut_pos:len(vcf_rec.alts[0])] != vcf_rec.alts[0]):
                        continue
                
                if ('S' in read.cigarstring and  # type: ignore
                    mean(read.query_alignment_qualities) < clip_qual_cutoff):  # type: ignore
                    continue
                
                if read.flag & 0x40:  # read first in pair
                    # ADD READ TO DICT OR SOMETHING
                    mut_reads[mut_sample_name].append(read)
                else:  # read second in pair
                    if read.flag & 0x10:
                        mate = bam_reader_dict[mut_sample_name].mate(read)
                        if mate.reference_end is None:
                            continue  # ?
                        if read.reference_start <= mate.reference_end:  # check with Peter, does this map on to his code accurately
                            read_start = mate.query_alignment_end + 1
                    else:
                        if read.reference_end >= read.next_reference_start:
                            read_end = read.next_reference_start - 1
                    if read_start <= vcf_rec.pos <= read_end:
                        # ADD READ TO DICT OR SOMETHING
                        mut_reads[mut_sample_name].append(read)
            ### end
        ### remove_dups_with_wobble
        for _, reads in mut_reads.items():
            if len(reads) == 0:
                continue
                # want the start of the record, the end
            ### start_mate_end_pairs()
            # incidentally, I suppose hairpin only works for paired data?
            sorted_ends = []
            for read in reads:
                mate = bam.mate(read)
                
                if any(x is None for x in [read.reference_start,
                                        read.reference_end,
                                        read.cigartuples,
                                        mate.reference_start,
                                        mate.reference_end,
                                        mate.cigartuples]):
                    continue  # ?
                
                # this gets pos wrt to reference, not query sequence, is that desired?
                start: int = read.reference_start
                end: int = read.reference_end
                mate_start: int  = mate.reference_start
                mate_end: int  = mate.reference_end
                
                # Peter
                # behaviour on cig = None?     
                # Julia XAM gets cigar info differently, checking CG:B,I tag
                # does this matter?
                if read.cigartuples[0][0] == Ops.soft:
                    start -= read.cigartuples[0][1]
                if read.cigartuples[-1][0] == Ops.soft:
                    end += cig[-1][1]
                if mate.cigartuples[0][0] == Ops.soft:
                    mate_start -= mate.cigartuples[0][1]
                if mate.cigartuples[-1][0] == Ops.soft:
                    mate_end += mate.cigartuples[-1][1]
                
                # appears mate posns simply aren't assigned if none
                sorted_ends.append(sorted([start, end, mate_start, mate_end]))
            ### end
            # I don't really understand this
            sorted_ends: list[list[int]] = sorted(sorted_ends)
            min_ends: list[list[int]] = [sorted_ends.pop(0)]
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
            del(i)
        ### end
        
        ### check hairpin filter
        mut_read_pos_f: list[int] = []
        mut_read_pos_r: list[int] = []
        mut_read_fracs_f: list[float] = []
        mut_read_fracs_r: list[float] = []
        aln_scores: list[float] = []
        for _, reads in mut_reads.items():
            for read in reads:
                if any([x is None for x in [read.reference_start, read.reference_end]]):
                    continue  # ?
                
                mut_pos = ref2seq.ref2querypos(read, vcf_rec.pos)
                if mut_pos == -1:
                    continue  # ?
                if read.flag & 0x10:
                    read_loc = read.reference_end - mut_pos + 1
                    mut_read_fracs_r.append(read_loc / (read.reference_start - read.reference_end + 1))
                    mut_read_pos_r.append(read_loc)
                else:
                    read_loc = (mut_pos - read.reference_start + 1)
                    mut_read_fracs_f.append(read_loc / (read.reference_end - read.reference_start + 1))
                    mut_read_pos_f.append(read_loc)
                try:
                    read.get_tag('AS')
                except KeyError:
                    continue  # ?
                aln_scores.append(read.get_tag('AS') / read.query_length)  # or should this be .query_alignment_length? (Peter)
        al_filt = median(aln_scores) <= al_thresh
        
        if len(mut_read_pos_f) > 1:
            mad_f = max(mut_read_pos_f) - min(mut_read_pos_f)
            sd_f = stdev(mut_read_pos_f)
        if len(mut_read_pos_r) > 1:
            mad_r = max(mut_read_pos_r) - min(mut_read_pos_r)
            sd_r = stdev(mut_read_pos_r)
        # hairpin conditions from Ellis et al.
        hp_filt = True
        # these branches all lead to the same result!
        if len(mut_read_pos_f) > 1 and len(mut_read_pos_r) > 1:
            frac_lt_thresh = sum([x <= cent90_thresh for x in mut_read_fracs_f + mut_read_fracs_r]) / (len(mut_read_pos_f) + len(mut_read_pos_r))
            if (frac_lt_thresh < 0.9 or
                (mad_f > 2 and mad_r > 2 and sd_f > 2 and sd_r > 2) or
                (mad_f > 1 and sd_f > 10) or
                (mad_r > 1 and sd_r > 10)):
                hp_filt = False
        elif len(mut_read_pos_f) > 1:
            if (((sum([x <= cent90_thresh for x in mut_read_pos_f]) / len(mut_read_pos_f)) < 0.9) and
                mad_f > 0 and
                sd_f > 4):
                hp_filt = False
        elif len(mut_read_pos_r) > 1:
            if (((sum([x <= cent90_thresh for x in mut_read_pos_r]) / len(mut_read_pos_r)) < 0.9) and
                mad_r > 0 and
                sd_r > 4):
                hp_filt = False
        else:
            hp_filt = False
        ### end 

        ### update vcf record
        if al_filt:
            vcf_rec.filter.add("ALF")
        if hp_filt:
            vcf_rec.filter.add("HPF")
        
        # try except    
        vcf_out.write(vcf_rec)

    return


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(prog="hairpin")
    parser.add_argument('-i', '--vcf-in', help="path to input vcf", nargs=1, type=str, required=True)
    parser.add_argument('-o', '--vcf-out', help="path to vcf out", nargs=1, type=str, required=True)
    parser.add_argument('-b', '--bams', help="list of paths to bams for samples in input vcf, whitespace separated", nargs='+', type=list, required=True)
    parser.add_argument('-cq', '--clip-quality-cutoff', default=35)
    parser.add_argument('-mq', '--min-mapping-quality', default=11)
    parser.add_argument('-mb', '--min-base-quality', default=25)
    parser.add_argument('-ms', '--max-read-span', default=6)
    parser.add_argument('-al', '--AL-filter-threshold', default=0.93)
    parser.add_argument('-c9', '--cent90-threshold', default=0.15)
    
    args = parser.parse_args()
    
    