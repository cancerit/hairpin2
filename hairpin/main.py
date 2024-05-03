import pysam
from enum import Enum

SOFT_CLIP_OP = 4
# for start_end_mate_pairs
Pairs = Enum('Pairs', [('start', 0), ('end', 1), ('mate_start', 2), ('mate_end', 3)])

# Peter modifies a dict, nonref_reads, inplace
# i.e. each sample has a dict entry
# dict is sample name : bam record
# return bam record, handle dict above
def get_mutant_reads(
    vcf_record, # called per record
    bam, # type depends on streaming approach
    min_basequal: int,
    clip_qual_cutoff: int,
    min_mapqual: int
): # return type probably pysam.AlignedSegment
    return

# pysam AlignmentFile.Fetch will return iterator over reads
# which yields AlignedSegment
def start_end_mate_pairs(
    record: pysam.AlignedSegment,
    bam: pysam.AlignmentFile
) -> list[int]:
    # want the start of the record, the end
    # start of mate, end of mate
    # .reference_start
    mate: pysam.AlignedSegment = bam.mate(record)
    # Julia XAM gets cigar info differently, checking CG:B,I tag
    # does this matter?
    cig: list[tuple[int, int]] = record.cigartuples
    mate_cig: list[tuple[int, int]] = mate.cigartuples
    
    start: int = record.reference_start
    end: int = record.reference_end
    mate_start: int = record.reference_start
    mate_end: int = record.reference_end
    
    # behaviour on cig = None?
    if cig[0][0] == SOFT_CLIP_OP:
        start -= cig[0][1]
    if cig[-1][0] == SOFT_CLIP_OP:
        end += cig[-1][1]
    if mate_cig[0][0] == SOFT_CLIP_OP:
        mate_start -= mate_cig[0][1]
    if mate_cig[-1][0] == SOFT_CLIP_OP:
        mate_end += mate_cig[-1][1]
    
    # appears mate posns simply aren't assigned if none
    return [start, end, mate_start, mate_end]


# Peter does this in place. Any particular reason?
# is structure of nonref_reads optimal/appropriate
def remove_dups_with_wobble(
    nonref_reads: dict,
    max_span_ends: int
) -> dict:
    
    return 'hello world'
    

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
    
    return

if __name__ == '__main__':
    # do stuff
    print('hello world')