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
import copy
import pytest


# BASIS PATH TESTING
# test every node and edge at least once
# N.B.
# pysam guards against:
# quality and seq length mismatch
# reference id is none
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
def test_path_unsupported_alt():
    with pytest.raises(ValueError):
        hp2.validate_read(read=r,
                          vcf_start=99,
                          vcf_stop=100,
                          alt='8',
                          mut_type='S',
                          min_mapqual=11,
                          min_clipqual=35,
                          min_basequal=25)


@pytest.mark.validate
def test_path_unsupported_mut_type():
    with pytest.raises(ValueError):
        hp2.validate_read(read=r,
                          vcf_start=99,
                          vcf_stop=100,
                          alt='A',
                          mut_type='8',
                          min_mapqual=11,
                          min_clipqual=35,
                          min_basequal=25)


@pytest.mark.validate
def test_path_missing_mc():
    rc = copy.deepcopy(r)
    rc.set_tag('MC', None)
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.READ_FIELDS_MISSING.value


@pytest.mark.validate
def test_path_missing_field():
    rc = copy.deepcopy(r)
    rc.cigarstring = None
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.READ_FIELDS_MISSING.value


@pytest.mark.validate
def test_path_set_flag_mapqual_clipqual():
    rc = copy.deepcopy(r)
    rc.flag = 0x200
    rc.cigarstring = '1S9M'
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=30,
                             min_clipqual=40,
                             min_basequal=25) == (c.ValidatorFlags.FLAG.value | c.ValidatorFlags.MAPQUAL.value | c.ValidatorFlags.CLIPQUAL.value)


@pytest.mark.validate
def test_path_sub_not_aligned():
    assert hp2.validate_read(read=r,
                             vcf_start=200,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.NOT_ALIGNED.value


@pytest.mark.validate
def test_path_bad_sub():
    assert hp2.validate_read(read=r,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='T',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=50) == (c.ValidatorFlags.NOT_ALT.value | c.ValidatorFlags.BASEQUAL.value)


@pytest.mark.validate
def test_path_good_sub():
    assert hp2.validate_read(read=r,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.CLEAR.value


# checks cigar ops
@pytest.mark.validate
def test_path_del_bad_op():
    assert hp2.validate_read(read=r,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='.',
                             mut_type='D',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.BAD_OP.value


# 2bp del
@pytest.mark.validate
def test_path_good_del():
    rc = copy.deepcopy(r)
    rc.cigarstring = '4M2D6M'
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=101,
                             alt='.',
                             mut_type='D',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.CLEAR.value


@pytest.mark.validate
def test_path_ins_not_aligned():
    assert hp2.validate_read(read=r,
                             vcf_start=200,
                             vcf_stop=100,
                             alt='A',
                             mut_type='I',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.NOT_ALIGNED.value


@pytest.mark.validate
def test_path_ins_short():
    assert hp2.validate_read(read=r,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='ATTTTTTTTTTTTTT',
                             mut_type='I',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.SHORT.value


@pytest.mark.validate
def test_path_bad_ins():
    assert hp2.validate_read(read=r,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='AC',
                             mut_type='I',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == (c.ValidatorFlags.BAD_OP.value | c.ValidatorFlags.NOT_ALT.value)


@pytest.mark.validate
def test_path_good_ins():
    rc = copy.deepcopy(r)
    rc.cigarstring = '5M2I3M'
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='AA',
                             mut_type='I',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.CLEAR.value


@pytest.mark.validate
def test_path_overlap():
    rc = copy.deepcopy(r)
    rc.flag = 0x83
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.OVERLAP.value


@pytest.mark.validate
def test_path_no_overlap():
    rc = copy.deepcopy(r)
    rc.flag = 0x83
    rc.set_tag('MC', '3M')
    assert hp2.validate_read(read=rc,
                             vcf_start=99,
                             vcf_stop=100,
                             alt='A',
                             mut_type='S',
                             min_mapqual=11,
                             min_clipqual=35,
                             min_basequal=25) == c.ValidatorFlags.CLEAR.value
