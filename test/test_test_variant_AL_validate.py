from hairpin2.main import test_variant_AL
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
    expected = c.ALFilter(flag=True, code=2, avg_as=0.5)
    s1r1 = copy.deepcopy(r)  # no AS, cover except KeyError
    s1r2 = copy.deepcopy(r)
    s1r2.set_tag('AS', 50)  # low AS
    result = test_variant_AL([s1r1, s1r2])
    assert expected == result


@pytest.mark.validate
def test_path_AL_false_code_2():
    expected = c.ALFilter(flag=False, code=2, avg_as=0.99)
    s1r1 = copy.deepcopy(r)
    s1r1.set_tag('AS', 99)  # high AS
    result = test_variant_AL([s1r1])
    assert expected == result


@pytest.mark.validate
def test_path_AL_false_code_3():
    expected = c.ALFilter(code=3)
    result = test_variant_AL([])
    assert expected == result
