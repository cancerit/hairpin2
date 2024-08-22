from hairpin2 import main as hp2
from hairpin2 import constants as c
import pysam

r = pysam.AlignedSegment()
r.query_name = 'read1'
r.query_sequence = 'CTGDAAAACC'
r.query_qualities = pysam.qualitystring_to_array('KKKKKKKKKK')
r.flag = 0x2
r.reference_id = 0
r.reference_start = 95
r.next_reference_start = 105
r.mapping_quality = 20
r.cigarstring = '10M'
r.cigartuples = [(0,10)]
r.set_tag('MC', '10M')


def test_validate_read_ideal():
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.CLEAR.value


def test_validate_read_missing_cigar():
    rc = r
    rc.cigartuples = None
    assert hp2.validate_read(read=r, vcf_start=99, vcf_stop=100, vcf_rlen=1, alt='A', min_mapqual=11, min_clipqual=35, min_basequal=25) == c.ValidatorFlags.READ_FIELDS_MISSING.value

