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


from hairpin2 import main as hp2
from hairpin2 import constants as c
import pysam
import pytest
import copy
from functools import partial


# BASIS PATH TESTING
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


# max spans...
### 21/10/24 HERE: I separated out max spans to make testing easier
### or at least to make understanding max spans easier
### I'm not super keen on the separate function so
### maybe it's going back in but now I actually understand it
### perhaps do test suite without it, put it back and cover at the end.
# don't forget to install updated main.py


# use this test to test all initial loops, and other tests to test unique outcomes
# e.g. read_flag != clear and so on
@pytest.mark.validate
def test_path_insufficient_reads():
    expected = c.Filters(AL=c.ALFilter(code=3),
                         HP=c.HPFilter(code=3))
    readd = {'S1': []}
    # using defaults where not otherwise noted
    actual = hp2.test_variant(
        vstart=166,
        vstop=167,
        alt='A',
        region_reads_by_sample=readd,
        mut_type='S',
        al_thresh=0.93,
        max_span=6,
        position_fraction_thresh=0.15,
        read_validator=partial(hp2.validate_read,
                               min_mapqual=11,
                               min_clipqual=35,
                               min_basequal=25)
    )
    assert expected == actual


# N.B. copy r to several reads
# combine into dict[str, Iterable[pysam.AlignedSegment]]
# where keys are samples
# yank this to register
@pytest.mark.validate
def test_path_xxx():
    pass


@pytest.mark.validate
def test_path_AL_true_code_2():
    expected = c.ALFilter(flag=True, code=2, avg_as=0.5)
    s1r1 = copy.deepcopy(r)  # no AS, cover except KeyError
    s1r2 = copy.deepcopy(r)
    s1r2.set_tag('AS', 50)  # low AS
    readd = {'S1': [s1r1, s1r2]}
    result = hp2.test_variant(
        vstart=166,
        vstop=167,
        alt='A',
        region_reads_by_sample=readd,
        mut_type='S',
        al_thresh=0.93,
        max_span=-1,  # don't trigger PCR dedup
        position_fraction_thresh=0.15,
        read_validator=partial(hp2.validate_read,
                               min_mapqual=11,
                               min_clipqual=35,
                               min_basequal=25)
    )
    assert expected == result.AL


@pytest.mark.validate
def test_path_AL_false_code_2_HP_false_code_3():
    expected = c.Filters(c.ALFilter(flag=False, code=2, avg_as=0.99),
                         c.HPFilter(code=3))
    s1r1 = copy.deepcopy(r)
    s1r1.set_tag('AS', 99)  # high AS
    readd = {'S1': [s1r1]}
    result = hp2.test_variant(
        vstart=166,
        vstop=167,
        alt='A',
        region_reads_by_sample=readd,
        mut_type='S',
        al_thresh=0.93,
        max_span=-1,  # don't trigger PCR dedup
        position_fraction_thresh=0.15,
        read_validator=partial(hp2.validate_read,
                               min_mapqual=11,
                               min_clipqual=35,
                               min_basequal=25)
    )
    assert expected == result


@pytest.mark.validate
def test_path_AL_false_code_3():
    expected = c.ALFilter(code=3)
    s1r1 = copy.deepcopy(r)
    readd = {'S1': [s1r1]}
    result = hp2.test_variant(
        vstart=166,
        vstop=167,
        alt='A',
        region_reads_by_sample=readd,
        mut_type='S',
        al_thresh=0.93,
        max_span=-1,  # don't trigger PCR dedup
        position_fraction_thresh=0.15,
        read_validator=partial(hp2.validate_read,
                               min_mapqual=11,
                               min_clipqual=35,
                               min_basequal=25)
    )
    assert expected == result.AL


@pytest.mark.validate
def test_path_HP_insufficient_reads():
    s1r1 = copy.deepcopy(r)
