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
    """
    Find query position on read for a given reference position.
    """
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
    if type(cig_str) != str:
        raise ValueError("cig_str not of type str")
    if type(ref_start) != int:
        raise ValueError("ref_start not of type int")
    if not cig_str:
        raise ValueError("cigar string empty")
    opl: list[tuple[int, str]] = []
    digit_accumulator: str = ""
    for char in cig_str:
        if char.isdigit():
            digit_accumulator += char
        elif char in ["M", "D", "N", "=", "X", "I", "P", "H", "S"]:
            opl.append((int(digit_accumulator), char))
            digit_accumulator = ""
        else:
            raise CigarError("could not interpret cigar string {}".format(cig_str))
    if not opl:
        raise CigarError("could not interpret cigar string {}".format(cig_str))
    for op_len, op_code in opl:
        if op_code in ["M", "D", "N", "=", "X"]:  # if op consumes ref
            ref_start += op_len
    return ref_start
