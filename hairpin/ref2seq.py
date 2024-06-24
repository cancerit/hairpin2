import pysam

def ref2querypos(
            record: pysam.AlignedSegment,
            ref_pos: int
) -> int | None:
    pos_aln = record.get_aligned_pairs()
    for aln_pair in pos_aln:
        if aln_pair[1] == ref_pos:
            return aln_pair[0]
    raise IndexError('reference position not covered by read')


def pos2op(
    seq_pos: int,
    record: pysam.AlignedSegment
) -> int:
    cig = record.cigartuples
    if cig is None:
        raise ValueError('no cigar tuples available for pysam record') # No cigar tuples for record
    sum_len = 0
    while True:
        cig_pair = cig.pop(0)
        sum_len += cig_pair[1]
        if seq_pos < sum_len:
            return cig_pair[0]

