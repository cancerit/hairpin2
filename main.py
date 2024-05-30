import pysam
from enum import Enum
from typing import Optional
from hairpin import ref2seq
import statistics

Ops = Enum('Ops',
           ['match', 'ins', 'delete', 'skip', 'soft', 'hard', 'pad', 'equal', 'diff', 'back'],
           start = 0)
# for start_end_mate_pairs
Pairs = Enum('Pairs',
             ['start', 'end', 'mate_start', 'mate_end'],
             star = 0)


# Peter does this in place. Any particular reason?
# is structure of nonref_reads optimal/appropriate
def remove_dups_with_wobble(
    nonref_reads: dict,
    max_span_ends: int
) -> dict:
    
    return dict()
    

# filters analysed on a cohort basis
def test_filters(
    vcf_posn: int,
    mutant_reads: list[pysam.AlignedSegment],
    cent90_thresh: float,
    AL_filt_thresh: float
) -> tuple[bool, bool]:
    
    mut_pos_f: list[int] = []
    mut_fracs_f: list[float] = []
    mut_pos_r: list[int] = []
    mut_fracs_r: list[float] = []
    aln_scores: list[float] = []
     
    return tuple()


# is streaming approach necessary?
def main(
    bam_paths: list,
    intervals, # type?
    vcf_in_path: str,
    vcf_out_path: str,
    clip_qual_cutoff: int,
    min_mapqual: int,
    min_basequal: int,
    max_span: int,
    AL_thresh: float,
    cent90_thresh:float,
    header: bool
) -> None:
    
    vcf_obj: pysam.VariantFile = pysam.VariantFile(vcf_in_path)
    sample_names: list[str] = list(vcf_obj.header.samples)
    mut_reads: dict[str, list[pysam.AlignedSegment]] = {key: [] for key in sample_names}
    
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
    
    # init output
    if header:
        pass
    
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
            for read in bam_reader_dict[mut_sample_name].fetch(vcf_rec.chrom, vcf_rec.start, vcf_rec.stop) # type: ignore
                
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
                    statistics.mean(read.query_alignment_qualities) < clip_qual_cutoff):  # type: ignore
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
        vcf_rec
        mut_reads
        samples_w_mutants
        
        for samp, reads in mut_reads:
            for read in reads:
                
    
    return

if __name__ == '__main__':
    # do stuff
    print('hello world')