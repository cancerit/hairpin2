import pysam
import pytest

from hairpin2.infrastructure.structures import ExtendedRead


@pytest.fixture(scope="session")
def read():
    """
    Creates a pytest fixture that generates and initialises a `pysam.AlignedSegment` object.

    This fixture provides a pre-configured aligned read object that mimics a sequencing
    read, initialises with attributes commonly used in genomic analyses. The fixture
    operates with a session scope, meaning it will be available and initialises once per
    test session.

    Fixture scope:
        session

    Returns:
        pysam.AlignedSegment: An aligned segment object initialises with predefined
        attributes.

    Raises:
        None
    """
    read = pysam.AlignedSegment()
    read.query_name = "read1"
    read.query_sequence = "CTGDAAAACC" * 10
    read.query_qualities = pysam.qualitystring_to_array("AAAAAAAAAA" * 10)
    read.flag = 0x43
    read.reference_id = 0
    read.reference_start = 100
    read.next_reference_start = 100
    read.mapping_quality = 20
    read.cigarstring = "100M"
    read.set_tag("MC", "100M")
    extr = ExtendedRead(read)
    return extr