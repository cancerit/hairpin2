import pysam
from enum import Enum

SOFT_CLIP_OP = 4
# for start_end_mate_pairs
Pairs = Enum('Pairs', [('start', 0), ('end', 1), ('mate_start', 2), ('mate_end', 3)])



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
    
    
        
    
        
    
    
    
    
    
    
    
    return [] 

if __name__ == '__main__':
    # do stuff
    print('hello world')