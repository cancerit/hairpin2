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
from copy import deepcopy
import pysam

from hairpin2.filters import DVF


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


### TEST NODES AND EDGES
def test_path_insufficient_reads():
    dv = DVF.Filter(fixed_params=DVF.Params())
    readsin: dict[str, list[pysam.AlignedSegment]] = {'_': []}
    readsout, result = dv.test('_', readsin)
    assert result.code == DVF.DVCodes.INSUFFICIENT_READS
    assert result.flag == None
    assert readsin == readsout


def test_path_duplicated():
    boundary_wobble = 10
    dv = DVF.Filter(fixed_params=DVF.Params((boundary_wobble - 1)))
    r1 = deepcopy(r)
    r1.reference_start -= boundary_wobble
    r1.next_reference_start -= boundary_wobble
    readsin: dict[str, list[pysam.AlignedSegment]] = {'_': [r, r, r1], '_2': [r]}
    readsout, result = dv.test('_', readsin)
    assert result.code == DVF.DVCodes.DUPLICATION
    assert result.flag == True
    assert readsin != readsout
    assert len(readsout.keys()) == 1
    assert len(readsout.values()) == 1


def test_path_duplicated_below_threshold():
    dv = DVF.Filter(fixed_params=DVF.Params(read_number_difference_threshold=1, nsamples_threshold=1))
    readsin: dict[str, list[pysam.AlignedSegment]] = {'_': [r, r]}
    readsout, result = dv.test('_', readsin)
    assert result.code == DVF.DVCodes.DUPLICATION
    assert result.flag == False
    assert readsin != readsout
    assert len(readsout.keys()) == 1
    assert len(readsout.values()) == 1


# test parameters behave as intended
def test_check_window():
    """
    duplication_window_size must define an inclusive window size in number of bases within which reads are considered duplicated
    """
    boundary_wobble = 10
    r1 = deepcopy(r)
    r1.reference_start -= boundary_wobble
    r1.next_reference_start -= boundary_wobble
    readsin: dict[str, list[pysam.AlignedSegment]] = {'_': [r, r1]}
    dv = DVF.Filter(fixed_params=DVF.Params((boundary_wobble - 1)))
    _, result = dv.test('_', readsin)
    assert result.flag == False
    dv = DVF.Filter(fixed_params=DVF.Params((boundary_wobble)))
    _, result = dv.test('_', readsin)
    assert result.flag == True


