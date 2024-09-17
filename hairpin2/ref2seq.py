# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import pysam
from hairpin2 import constants as c

def ref2querypos(
            bam_record: pysam.AlignedSegment,
            ref_pos: int,
            get_cig: bool = True
) -> tuple[int, int | None]:
    pos_aln = bam_record.get_aligned_pairs()
    query_pos = pos_op = None
    for aln_pair in pos_aln:
        if aln_pair[1] == ref_pos:
            query_pos = aln_pair[0]
    if query_pos is None or len(pos_aln) == 0:
        raise IndexError('reference position not covered by read')
    elif get_cig:
        dist2op = ref_pos - bam_record.reference_start + 1  # since position is 0-indexed, add 1 to get distance
        cig = bam_record.cigartuples
        if cig is None or len(cig) == 0:
            raise ValueError('no cigar tuples available for pysam record')
        sum_len = 0
        while len(cig) > 0:
            cig_pair = cig.pop(0)
            if cig_pair[0] != c.Ops.SOFT.value:
                sum_len += cig_pair[1]
                if dist2op <= sum_len:
                    pos_op = cig_pair[0]
        if pos_op is None:
            raise ValueError('cigar op could not be recovered')
    return query_pos, pos_op


def ref_end_via_cigar(
    cig_str: str,
    ref_start: int
) -> int:
    if not cig_str[0].isdigit() or len(cig_str) < 2:
        raise ValueError('cigar string misformatted')
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
        if op_code in ['M','D','N','=','X']:
            ref_start += int(op_len)
    return ref_start
                
