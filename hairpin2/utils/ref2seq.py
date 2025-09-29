# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import pysam


def ref2querypos(
    read: pysam.AlignedSegment,
    ref_pos: int,
) -> int:
    pos_aln = read.get_aligned_pairs()
    query_pos = None
    for aln_pair in pos_aln:
        if aln_pair[1] == ref_pos:
            query_pos = aln_pair[0]
    if query_pos is None or len(pos_aln) == 0:
        raise ValueError("reference position not covered by read")
    return query_pos


class CigarError(ValueError):
    pass


def ref_end_via_cigar(cig_str: str, ref_start: int) -> int:
    if (
        not cig_str[0].isdigit()
        or not all([(c.isalnum() or c == "=") for c in cig_str])
        or len(cig_str) < 2
    ):
        raise CigarError("could not interpret cigar string {}".format(cig_str))
    cig_l = []
    digit_accumulator: str = ""
    for char in cig_str:
        if char.isdigit():
            digit_accumulator += char
        else:
            cig_l.append(digit_accumulator)
            cig_l.append(char)
            digit_accumulator = ""
    cig_t = list(zip(cig_l[0::2], cig_l[1::2]))
    for op_len, op_code in cig_t:
        if op_code in ["M", "D", "N", "=", "X"]:
            ref_start += int(op_len)
    return ref_start
