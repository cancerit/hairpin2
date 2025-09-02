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
from copy import deepcopy
from typing import cast
from pysam import AlignedSegment, VariantRecord, qualitystring_to_array
from hairpin2.flaggers import DVF
from hairpin2.structures import ReadView
from .helpers import comp_ReadView, unsafe_construct


# perfect read pair:
r = AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC' * 10
r.query_qualities = qualitystring_to_array('AAAAAAAAAA' * 10)
r.flag = 0x43
r.reference_id = 0
r.reference_start = 100
r.next_reference_start = 100
r.mapping_quality = 20
r.cigarstring = '100M'
r.set_tag('MC', '100M')


fake_rec = cast(VariantRecord, cast(object, ""))  # not used by DVF test method
fake_alt = "_"  # not material for DVF test record, just for record keeping


### TEST NODES AND EDGES
def test_path_insufficient_reads():
    dv = DVF.FlaggerDVF(fixed_params=DVF.FixedParamsDVF())
    reads = ReadView({'_': []})
    rsnapshot = deepcopy(reads)
    test_data = unsafe_construct(DVF.VarParamsDVF, record=None, alt=fake_alt, reads=reads)
    dv._var_params = cast(DVF.VarParamsDVF, test_data)  # unsafe prime
    
    result = dv.test()
    assert result.code == DVF.CodesDVF.INSUFFICIENT_READS
    assert result.flag == None
    assert comp_ReadView(reads, rsnapshot)


def test_path_duplicated():
    dv = DVF.FlaggerDVF(fixed_params=DVF.FixedParamsDVF())
    r1 = deepcopy(r)
    r1.is_duplicate = True
    reads = ReadView({'_': [r, r1, r1]}, "no_validate")
    rsnapshot = deepcopy(reads)
    test_data = unsafe_construct(DVF.VarParamsDVF, record=None, alt=fake_alt, reads=reads)
    dv._var_params = cast(DVF.VarParamsDVF, test_data)

    result = dv.test()
    assert result.code == DVF.CodesDVF.DUPLICATION
    assert result.flag == True
    assert comp_ReadView(reads, rsnapshot)


# TODO: test loss ratio, nsamples params

