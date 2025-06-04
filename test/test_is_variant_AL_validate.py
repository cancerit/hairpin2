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

from hairpin2.main import variant_AL_threshold
from hairpin2 import constants as c
import pysam
import pytest
import copy


# BASIS PATH TESTING (ish)
# test every node and edge at least once
# ----
# perfect read pair:
r = pysam.AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC' * 10
r.query_qualities = pysam.qualitystring_to_array('AAAAAAAAAA' * 10)
r.flag = 0x43
r.reference_id = 0
r.reference_start = 100
r.next_reference_start = 100
r.mapping_quality = 20
r.cigarstring = '100M'
r.set_tag('MC', '100M')


@pytest.mark.validate
def test_path_AL_true_code_2():
    al = c.ALFilter()
    expected = c.ALFilter(flag=True, code=c.ALFCodes.ON_THRESHOLD, avg_as=0.5)
    s1r1 = copy.deepcopy(r)  # no AS, cover except KeyError
    s1r2 = copy.deepcopy(r)
    s1r2.set_tag('AS', 50)  # low AS
    variant_AL_threshold(al, [s1r1, s1r2])
    assert expected == al


@pytest.mark.validate
def test_path_AL_false_code_2():
    al = c.ALFilter()
    expected = c.ALFilter(flag=False, code=c.ALFCodes.ON_THRESHOLD, avg_as=0.99)
    s1r1 = copy.deepcopy(r)
    s1r1.set_tag('AS', 99)  # high AS
    variant_AL_threshold(al, [s1r1])
    assert expected == al


@pytest.mark.validate
def test_path_AL_false_code_3():
    al = c.ALFilter()
    expected = c.ALFilter(flag=False, code=c.ALFCodes.INSUFFICIENT_SUPPORT)
    variant_AL_threshold(al, [])
    assert expected == al
