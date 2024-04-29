import pysam


# pysam AlignmentFile.Fetch will return iterator over reads
# which I think itself returns iterator row
# hence the typing of this function
# needs checking
def start_end_mate_pairs(record: pysam.IteratorRow) -> list[int]:
    
    return [] 

if __name__ == '__main__':
    