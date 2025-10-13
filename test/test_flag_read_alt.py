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
import copy

import pysam
from hairpin2.readqc import ValidatorFlags, qc_read_alt_specific
from pysam.libcalignedsegment import SAM_FLAGS as s
from pytest import raises

# perfect read pair:
r = pysam.AlignedSegment()
r.query_name = "read1"
r.query_sequence = "CTGDAAAACC"
r.query_qualities = pysam.qualitystring_to_array("AAAAAAAAAA")
r.flag = s.FPAIRED | s.FPROPER_PAIR | s.FREAD1  # 0x43
r.reference_id = 0
r.reference_start = 95
r.next_reference_start = 95
r.mapping_quality = 20
r.cigarstring = "10M"
r.set_tag("MC", "10M")


def test_path_unsupported_mut_type():
    with raises(ValueError):
        _ = qc_read_alt_specific(
            read=r,
            vcf_start=99,
            vcf_stop=100,
            alt="A",
            mut_type="8",  # pyright: ignore[reportArgumentType]
            min_basequal=25,
        )


def test_path_sub_not_aligned():
    expected = ValidatorFlags.NOT_ALIGNED
    result = qc_read_alt_specific(
        read=r, vcf_start=200, vcf_stop=100, alt="A", mut_type="S", min_basequal=25
    )
    assert expected == result


def test_path_bad_sub():
    expected = ValidatorFlags.NOT_ALT | ValidatorFlags.BASEQUAL
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=100, alt="T", mut_type="S", min_basequal=50
    )
    assert expected == result


def test_path_good_sub():
    expected = ValidatorFlags.CLEAR
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=100, alt="A", mut_type="S", min_basequal=25
    )
    assert expected == result


def test_path_del_short():
    expected = ValidatorFlags.SHORT
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=110, alt=".", mut_type="D", min_basequal=25
    )
    assert expected & result


def test_path_del_bad_op():
    expected = ValidatorFlags.BAD_OP
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=101, alt="C", mut_type="D", min_basequal=25
    )
    assert expected == result


# 2bp del
def test_path_good_del():
    expected = ValidatorFlags.CLEAR
    rc = copy.deepcopy(r)
    rc.cigarstring = "4M2D6M"
    result = qc_read_alt_specific(
        read=rc, vcf_start=98, vcf_stop=101, alt="CC", mut_type="D", min_basequal=25
    )
    assert expected == result


def test_path_ins_not_aligned():
    expected = ValidatorFlags.NOT_ALIGNED
    result = qc_read_alt_specific(
        read=r, vcf_start=200, vcf_stop=100, alt="A", mut_type="I", min_basequal=25
    )
    assert expected == result


def test_path_ins_short():
    expected = ValidatorFlags.SHORT
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=100, alt="ATTTTTTTTTTTTTT", mut_type="I", min_basequal=25
    )
    assert expected == result


def test_path_bad_ins():
    expected = ValidatorFlags.BAD_OP | ValidatorFlags.NOT_ALT
    result = qc_read_alt_specific(
        read=r, vcf_start=99, vcf_stop=100, alt="AC", mut_type="I", min_basequal=25
    )
    assert expected == result


def test_path_good_ins():
    expected = ValidatorFlags.CLEAR
    rc = copy.deepcopy(r)
    rc.cigarstring = "5M2I3M"
    result = qc_read_alt_specific(
        read=rc, vcf_start=99, vcf_stop=100, alt="AA", mut_type="I", min_basequal=25
    )
    assert expected == result
