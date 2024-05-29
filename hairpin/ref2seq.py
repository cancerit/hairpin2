import pysam

def ref2querypos(
            record: pysam.AlignedSegment,
            ref_pos: int
) -> int:
    pos_aln = record.get_aligned_pairs()
    while True:
        try:
            aln_pair = pos_aln.pop()
        except IndexError:
            return -1 # ref_pos not on read
        if aln_pair[1] == ref_pos:
            return aln_pair[0]


def pos2op(
    seq_pos: int,
    record: pysam.AlignedSegment
) -> int:
    cig = record.cigartuples
    if cig is None:
        exit(1)  # No cigar tuples for record
    sum_len = 0
    while True:
        try:
            cig_pair = cig.pop(0)
        except IndexError:
            raise RuntimeError  # seq_pos not in cigar string
        sum_len += cig_pair[1]
        if seq_pos < sum_len:
            return cig_pair[0]


"""
In Julia, ref2seq takes a run of cigar ops
and a position on the reference genome.
Here, the cigar ops come from the bam read/record
under examination, and the pos on ref comes from
the position given by the vcf for the vcf record
under examination (one or more bam records will be
examined for each vcf record). It returns an array
of 2 values, the position on the bam read, and the
cigar operation applicable to that position.

since peter wants to discard if not OP_MATCH
i.e. if clipped
we can use reference start in pysam
which is position where read begins alignment sans clipping etc
"""