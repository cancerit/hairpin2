# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import pysam


def ref2querypos(
            bam_record: pysam.AlignedSegment,
            ref_pos: int,
) -> int:
    pos_aln = bam_record.get_aligned_pairs()
    query_pos = None
    for aln_pair in pos_aln:
        if aln_pair[1] == ref_pos:
            query_pos = aln_pair[0]
    if query_pos is None or len(pos_aln) == 0:
        raise IndexError('reference position not covered by read')
    return query_pos


def ref_end_via_cigar(
    cig_str: str,
    ref_start: int
) -> int:
    if (not cig_str[0].isdigit() or
        not all([(c.isalnum() or c == '=') for c in cig_str]) or
            len(cig_str) < 2):
        raise ValueError('could not interpret cigar string {}'.format(cig_str))
    cig_l = []
    digit_accumulator: str = ''
    for char in cig_str:
        if char.isdigit():
            digit_accumulator += char
        else:
            cig_l.append(digit_accumulator)
            cig_l.append(char)
            digit_accumulator = ''
    cig_t = list(zip(cig_l[0::2], cig_l[1::2]))
    for op_len, op_code in cig_t:
        if op_code in ['M', 'D', 'N', '=', 'X']:
            ref_start += int(op_len)
    return ref_start
