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
from hairpin2.filters import ADF
import pysam
import copy


# TODO: parameter boundary testing (and elsewhere)
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


def test_path_insufficient_reads():
    ad = ADF.Filter(fixed_params=ADF.Params())
    readsin = []
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == None
    assert result.code == ADF.ADCodes.INSUFFICIENT_READS
    assert readsin == readsout


def test_path_insufficient_reads_2():
    ad = ADF.Filter(fixed_params=ADF.Params())
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    readsin = [r, rr]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == None
    assert result.code == ADF.ADCodes.INSUFFICIENT_READS
    assert readsin == readsout


def test_path_f_60ai_set():
    ad = ADF.Filter(fixed_params=ADF.Params())
    readsin =  [r, r]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == True
    assert result.code == ADF.ADCodes.SIXTYAI
    assert readsin == readsout


def test_path_f_60ai_noset():
    ad = ADF.Filter(fixed_params=ADF.Params())
    r1 = copy.deepcopy(r)
    r1.reference_start = 90
    readsin = [r, r1]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == False
    assert result.code == ADF.ADCodes.SIXTYAI
    assert readsin == readsout


def test_path_r_60ai_set():
    ad = ADF.Filter(fixed_params=ADF.Params())
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    readsin = [rr, rr]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == True
    assert result.code == ADF.ADCodes.SIXTYAI
    assert readsin == readsout


def test_path_r_60ai_noset():
    ad = ADF.Filter(fixed_params=ADF.Params())
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    rr1 = copy.deepcopy(rr)
    rr1.reference_start = 90
    readsin = [rr, rr1]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == False
    assert result.code == ADF.ADCodes.SIXTYAI
    assert readsin == readsout


def test_path_60bi_set():
    ad = ADF.Filter(fixed_params=ADF.Params())
    r1 = copy.deepcopy(r)
    r1.reference_start = 190
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    readsin = [r1, r1, rr, rr]
    readsout, result = ad.test('_', 198, readsin)
    assert result.flag == True
    assert result.code == ADF.ADCodes.SIXTYBI
    assert readsin == readsout


def test_path_60bi_noset():
    ad = ADF.Filter(fixed_params=ADF.Params())
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    readsin = [r, r, rr, rr]
    readsout, result = ad.test('_', 150, readsin)
    assert result.flag == False
    assert result.code == ADF.ADCodes.SIXTYBI
    assert readsin == readsout
