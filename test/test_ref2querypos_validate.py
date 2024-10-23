from hairpin2.ref2seq import ref2querypos
import pytest
import pysam


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
def test_path_indexerror():
    with pytest.raises(IndexError):
        ref2querypos(r, 1000)


@pytest.mark.validate
def test_path_good():
    expected = 5
    result = ref2querypos(r, 100)
    assert expected == result
