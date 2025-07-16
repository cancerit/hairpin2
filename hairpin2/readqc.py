# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
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
from pysam import AlignedSegment
from hairpin2 import ref2seq as r2s, constants as cnst
from statistics import mean
from typing import Literal


def qc_read_broad(
    read: AlignedSegment,
    vcf_start: int,
    min_mapqual: int,
    min_clipqual: int,
) -> cnst.ValidatorFlags:
    """
    When testing variants, test a read for various general - i.e. not specific to a
    particular variant - features identifying the read as a poor source of support
    for a variant. Used to disqualify reads for use in testing a variant.
    
    does not check for the following, as pysam already guards against:
        - quality and seq length mismatch
        - reference id is none
    """
    invalid_flag = cnst.ValidatorFlags.CLEAR  # 0 - evaluates false

    try:
        mate_cig = str(read.get_tag('MC'))
    except KeyError:
        mate_cig = None
    if any(flg is None for flg in
            [read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig]):
        invalid_flag |= cnst.ValidatorFlags.READ_FIELDS_MISSING
    else:
        if not (read.flag & 0x2) or read.flag & 0xE00:
            invalid_flag |= cnst.ValidatorFlags.FLAG

        if read.mapping_quality < min_mapqual:
            invalid_flag |= cnst.ValidatorFlags.MAPQUAL

        if ('S' in read.cigarstring and  # pyright: ignore[reportOperatorIssue]
                mean(read.query_alignment_qualities) < min_clipqual):  # pyright: ignore[reportUnknownMemberType, reportArgumentType]
            invalid_flag |= cnst.ValidatorFlags.CLIPQUAL

        # avoid analysing both read1 and mate if they both cover the variant
        if (not (invalid_flag & cnst.ValidatorFlags.FLAG)
                and not (read.flag & 0x40)):
            read_range = range(read.reference_start,
                               read.reference_end)  # pyright: ignore[reportArgumentType]
            mate_range = range(read.next_reference_start,
                               r2s.ref_end_via_cigar(mate_cig,  # pyright: ignore[reportArgumentType]
                                                     read.next_reference_start))
            ref_overlap = set(read_range).intersection(mate_range)
            if vcf_start in ref_overlap:
                invalid_flag |= cnst.ValidatorFlags.OVERLAP

    return invalid_flag


def qc_read_alt_specific(
    read: AlignedSegment,
    vcf_start: int,
    vcf_stop: int,
    alt: str,
    mut_type: Literal['S', 'D', 'I'],
    min_basequal: int
) -> cnst.ValidatorFlags:
    """
    When testing a variant, test a read for various features specific to the variant
    at hand which would identify that read as a poor source of support for that
    variant. Used to disqualify reads for use in testing a variant.
    """
    if mut_type not in ['S', 'D', 'I']:
        raise ValueError(
            'unsupported mut_type: {} - supports \'S\' (SUB) \'D\' (DEL) \'I\' (INS)'.format(mut_type))
    if read.query_sequence is None:
        raise ValueError(
            'read must have query sequence'
        )
    if read.query_qualities is None:
        raise ValueError(
            'read must have query qualities'
        )

    invalid_flag = cnst.ValidatorFlags.CLEAR

    if mut_type in ['S', 'I']:
        try:
            mut_pos = r2s.ref2querypos(read, vcf_start)
        except ValueError:
            invalid_flag |= cnst.ValidatorFlags.NOT_ALIGNED
        else:
            if mut_type == 'S':  # SUB
                if read.query_sequence[mut_pos:mut_pos + len(alt)] != alt:
                    invalid_flag |= cnst.ValidatorFlags.NOT_ALT
                if any([bq < min_basequal
                        for bq
                        in read.query_qualities[mut_pos:mut_pos + len(alt)]]):
                    invalid_flag |= cnst.ValidatorFlags.BASEQUAL
            if mut_type == 'I':  # INS - mut_pos is position immediately before insertion
                if mut_pos + len(alt) > read.query_length:
                    invalid_flag |= cnst.ValidatorFlags.SHORT
                else:
                    mut_alns = [(q, r)
                                for q, r
                                in read.get_aligned_pairs()
                                if q in range(mut_pos + 1, mut_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        invalid_flag |= cnst.ValidatorFlags.BAD_OP
                    if read.query_sequence[mut_pos + 1:mut_pos + len(alt) + 1] != alt:
                        invalid_flag |= cnst.ValidatorFlags.NOT_ALT
    elif mut_type == 'D':  # DEL
        rng = list(range(vcf_start, vcf_stop + 1))
        mut_alns = [q
                    for q, r
                    in read.get_aligned_pairs()
                    if r in rng]
        if len(mut_alns) != len(rng):
            invalid_flag |= cnst.ValidatorFlags.SHORT
        if (any([x is not None for x in mut_alns[1:-1]]) or
            any([x is None for x in [mut_alns[0], mut_alns[-1]]])):
                invalid_flag |= cnst.ValidatorFlags.BAD_OP

    return invalid_flag
