from hairpin2 import main as hp2
from hairpin2 import constants as c
import pysam
import copy

# pysam guards against:
# quality and seq length mismatch
# flag not set
# reference id is none

r = pysam.AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC'
r.query_qualities = pysam.qualitystring_to_array('AAAAAAAAAA')
r.flag = 0x2
r.reference_id = 0
r.reference_start = 95
r.next_reference_start = 105
r.mapping_quality = 20
r.cigarstring = '10M'
r.set_tag('MC', '10M')


# ideally there would be a test for each time read_flag is set
# i.e. test every path of achieving a given flag
# so far there's a test for each flag

def test_ideal():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.CLEAR.value


def test_mapqual():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=30, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.MAPQUAL.value


def test_not_aligned():
    assert hp2.validate_read(read=r, vcf_start=200, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.NOT_ALIGNED.value


def test_not_alt():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='T', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.NOT_ALT.value


def test_basequal():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=40) == c.ValidatorFlags.BASEQUAL.value


def test_read_short():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='ATTTTTTTTTTTTTT', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.SHORT.value




def test_missing_cigar():
    rc = copy.deepcopy(r)
    rc.cigarstring = None
    rc.cigartuples = None
    assert hp2.validate_read(read=rc, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.READ_FIELDS_MISSING.value


def test_bad_op():
    rc = copy.deepcopy(r)
    rc.cigartuples = [(c.Ops.EQUAL.value, 10)]
    assert hp2.validate_read(read=rc, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.BAD_OP.value


def test_clipqual():
    rc = copy.deepcopy(r)
    rc.cigarstring = '1S9M'
    assert hp2.validate_read(read=rc, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=40, min_basequal=25) == c.ValidatorFlags.CLIPQUAL.value


def test_no_overlap_0x10():
    rc = copy.deepcopy(r)
    rc.flag |= 0x10
    rc.set_tag('MC', '3M')
    assert hp2.validate_read(read=rc, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.NO_OVERLAP.value


def test_no_overlap():
    rc = copy.deepcopy(r)
    rc.next_reference_start = 98
    assert hp2.validate_read(read=rc, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.NO_OVERLAP.value


