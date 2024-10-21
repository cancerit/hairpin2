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

# S1 needs 2 good reads to give True on if len(mut_reads...) > 1
s1r1 = copy.deepcopy(r)
s1r2 = copy.deepcopy(r)
s1r3 = copy.deepcopy(r)
# S2 needs a bad read to give a False on if read_flag == ... CLEAR
# and len(mut_reads...)
s2r1 = copy.deepcopy(r)
s2r2 = copy.deepcopy(r)
s2r2.flag = 0xE00
iter1 = [s1r1, s1r2]
iter2 = [s2r1, s2r2]
readd = {'S1': iter1, 'S2': iter2}


# max spans...
### 21/10/24 HERE: I separated out max spans to make testing easier
### or at least to make understanding max spans easier
### I'm not super keen on the separate function so
### maybe it's going back in but now I actually understand it
### perhaps do test suite without it, put it back and cover at the end.
# don't forget to install updated main.py
@pytest.mark.validate
def test_path_simple():
    f = hp2.test_variant(
        vstart=160,
        vstop=161,
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
    breakpoint()
