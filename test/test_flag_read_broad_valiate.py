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


from hairpin2.main import flag_read_broad
from hairpin2 import constants as c
import pysam
import copy
import pytest


# BASIS PATH TESTING (ish)
# test every node and edge at least once
# ----
# perfect read pair:
r = pysam.AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC'
r.query_qualities = pysam.qualitystring_to_array('AAAAAAAAAA')
r.flag = 0x43
r.reference_id = 0
r.reference_start = 95
r.next_reference_start = 95
r.mapping_quality = 20
r.cigarstring = '10M'
r.set_tag('MC', '10M')


@pytest.mark.validate
def test_path_clear():
    expected = c.ValidatorFlags.CLEAR.value
    result = flag_read_broad(read=r,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


@pytest.mark.validate
def test_path_missing_mc():
    expected = c.ValidatorFlags.READ_FIELDS_MISSING.value
    rc = copy.deepcopy(r)
    rc.set_tag('MC', None)
    result = flag_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


@pytest.mark.validate
def test_path_missing_field():
    expected = c.ValidatorFlags.READ_FIELDS_MISSING.value
    rc = copy.deepcopy(r)
    rc.cigarstring = None
    result = flag_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=35)
    assert expected == result


@pytest.mark.validate
def test_path_set_flag_mapqual_clipqual():
    expected = (c.ValidatorFlags.FLAG.value
                | c.ValidatorFlags.MAPQUAL.value
                | c.ValidatorFlags.CLIPQUAL.value)
    rc = copy.deepcopy(r)
    rc.flag = 0x200
    rc.cigarstring = '1S9M'
    result = flag_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=30,
                             min_clipqual=40)
    assert expected == result


@pytest.mark.validate
def test_path_overlap():
    expected = c.ValidatorFlags.OVERLAP.value
    rc = copy.deepcopy(r)
    rc.flag = 0x83
    result = flag_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=40)
    assert expected == result


@pytest.mark.validate
def test_path_no_overlap():
    expected = c.ValidatorFlags.CLEAR.value
    rc = copy.deepcopy(r)
    rc.flag = 0x83
    rc.set_tag('MC', '3M')
    result = flag_read_broad(read=rc,
                             vcf_start=99,
                             min_mapqual=11,
                             min_clipqual=40)
    assert expected == result
