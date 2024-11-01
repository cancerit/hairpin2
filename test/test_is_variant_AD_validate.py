from hairpin2.main import is_variant_AD
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
def test_path_insufficient_reads():
    expected = c.ADFilter(code=3)
    result = is_variant_AD(0, [])
    assert expected == result


@pytest.mark.validate
def test_path_f_60ai_set():
    expected = c.ADFilter(flag=True, code=0)
    result = is_variant_AD(150, [r, r])
    assert expected == result


@pytest.mark.validate
def test_path_f_60ai_noset():
    expected = c.ADFilter(code=0)
    r1 = copy.deepcopy(r)
    r1.reference_start = 90
    result = is_variant_AD(150, [r, r1])
    assert expected == result


@pytest.mark.validate
def test_path_r_60ai_set():
    expected = c.ADFilter(flag=True, code=0)
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    result = is_variant_AD(150, [rr, rr])
    assert expected == result


@pytest.mark.validate
def test_path_r_60ai_noset():
    expected = c.ADFilter(code=0)
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    rr1 = copy.deepcopy(rr)
    rr1.reference_start = 90
    result = is_variant_AD(150, [rr, rr1])
    assert expected == result


@pytest.mark.validate
def test_path_60bi_set():
    expected = c.ADFilter(flag=True, code=1)
    r1 = copy.deepcopy(r)
    r1.reference_start = 190
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    result = is_variant_AD(198, [r1, r1, rr, rr])
    assert expected == result


@pytest.mark.validate
def test_path_60bi_noset():
    expected = c.ADFilter(code=1)
    rr = copy.deepcopy(r)
    rr.flag = rr.flag | 0x10
    result = is_variant_AD(150, [r, r, rr, rr])
    assert expected == result
