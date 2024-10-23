from hairpin2.ref2seq import ref_end_via_cigar
import pytest


# BASIS PATH TESTING (ish)
# test every node and edge at least once
# ----
@pytest.mark.validate
def test_path_valueerror():
    with pytest.raises(ValueError):
        ref_end_via_cigar('30M2I* ,', 0)


@pytest.mark.validate
def test_path_normal():
    expected = 50
    result = ref_end_via_cigar('20M10I30M', 0)
    assert expected == result
