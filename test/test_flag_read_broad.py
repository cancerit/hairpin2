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
from hairpin2.readqc import qc_read_broad, ValidatorFlags
import pysam
from pysam.libcalignedsegment import SAM_FLAGS as s
import copy


# perfect read pair:
r = pysam.AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC'
r.query_qualities = pysam.qualitystring_to_array('AAAAAAAAAA')
r.flag = s.FPAIRED | s.FPROPER_PAIR | s.FREAD1  # 0x43
r.reference_id = 0
r.reference_start = 95
r.next_reference_start = 95
r.mapping_quality = 20
r.cigarstring = '10M'
r.set_tag('MC', '10M')


def test_path_clear():
    expected = ValidatorFlags.CLEAR
    result = qc_read_broad(read=r,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


def test_path_missing_mc():
    expected = ValidatorFlags.READ_FIELDS_MISSING
    rc = copy.deepcopy(r)
    rc.set_tag('MC', None)
    result = qc_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


def test_path_missing_field():
    expected = ValidatorFlags.READ_FIELDS_MISSING
    rc = copy.deepcopy(r)
    rc.cigarstring = None
    result = qc_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


def test_path_set_flag_mapqual_clipqual():
    expected = (ValidatorFlags.FLAG
                | ValidatorFlags.MAPQUAL
                | ValidatorFlags.CLIPQUAL)
    rc = copy.deepcopy(r)
    rc.flag = s.FQCFAIL  # 0x200
    rc.cigarstring = '1S9M'
    result = qc_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=30,
                             min_clipqual=40)
    assert expected == result


def test_path_overlap():
    expected = ValidatorFlags.OVERLAP
    rc = copy.deepcopy(r)
    rc.flag = s.FPAIRED | s.FPROPER_PAIR | s.FREAD2  # 0x80
    result = qc_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=40)
    assert expected == result


def test_path_no_overlap():
    expected = ValidatorFlags.CLEAR
    rc = copy.deepcopy(r)
    rc.flag = s.FPAIRED | s.FPROPER_PAIR | s.FREAD2  # 0x80
    rc.set_tag('MC', '3M')
    result = qc_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=40)
    assert expected == result
