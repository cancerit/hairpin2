# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from copy import deepcopy

# TODO: user should only need to import one thing!
from hairpin2.read_preprocessors.dupmark import (
    DupmarkPreprocessor,
    FixedParamsDupmark,
    PrefilterParamsDupmark,
    RunParamsDupmark,
)
from hairpin2.structures import ReadView
from pysam import AlignedSegment, qualitystring_to_array

from test.helpers import unsafe_construct_params

# perfect read pair:
r = AlignedSegment()
r.query_name = "read1"
r.query_sequence = "CTGDAAAACC" * 10
r.query_qualities = qualitystring_to_array("AAAAAAAAAA" * 10)
r.flag = 0x43
r.reference_id = 0
r.reference_start = 100
r.next_reference_start = 100
r.mapping_quality = 20
r.cigarstring = "100M"
r.set_tag("MC", "100M")


def test_path_insufficient_reads():
    reads = ReadView({"_": []})
    assert sum(rd.has_tag("zD") for rd in reads.all) == 0

    # TODO: make default window size more available
    dp = DupmarkPreprocessor(
        prefilter_params=PrefilterParamsDupmark(min_mapq=0, min_avg_clipq=0, min_baseq=0),
        fixed_params=FixedParamsDupmark(duplication_window_size=6),
    )
    dp._var_params = unsafe_construct_params(RunParamsDupmark, reads=reads)

    dp.process()

    assert sum(rd.has_tag("zD") for rd in reads.all) == 0


def test_path_duplication():
    boundary_wobble = 10
    rdup = deepcopy(r)
    rnondup = deepcopy(r)  # create a non-dup
    rnondup.reference_start -= boundary_wobble
    rnondup.next_reference_start -= boundary_wobble
    # reads = [r, rdup, rnondup]

    reads = ReadView({"_": [r, rdup, rnondup]})
    assert sum(rd.has_tag("zD") for rd in reads.all) == 0

    dp = DupmarkPreprocessor(
        prefilter_params=PrefilterParamsDupmark(min_mapq=0, min_avg_clipq=0, min_baseq=0),
        fixed_params=FixedParamsDupmark(duplication_window_size=boundary_wobble - 1),
    )
    dp._var_params = unsafe_construct_params(RunParamsDupmark, reads=reads)

    dp.process()

    # _mark_stutter_dups(reads, boundary_wobble-1)  # TODO: check id in here!
    assert sum(rd.has_tag("zD") for rd in reads.all) == 1


# test parameters behave as intended
def test_check_window():
    """
    duplication_window_size must define an inclusive window size in number of bases within which reads are considered duplicated
    """
    boundary_wobble = 10
    r1 = deepcopy(r)
    r1.reference_start -= boundary_wobble
    r1.next_reference_start -= boundary_wobble

    reads = ReadView({"_": [deepcopy(r), deepcopy(r1)]})
    assert sum(rd.has_tag("zD") for rd in reads.all) == 0

    dp = DupmarkPreprocessor(
        prefilter_params=PrefilterParamsDupmark(min_mapq=0, min_avg_clipq=0, min_baseq=0),
        fixed_params=FixedParamsDupmark(duplication_window_size=boundary_wobble - 1),
    )
    dp._var_params = unsafe_construct_params(RunParamsDupmark, reads=reads)

    dp.process()
    assert sum(rd.has_tag("zD") for rd in reads.all) == 0

    reads = ReadView({"_": [deepcopy(r), deepcopy(r1)]})
    assert sum(rd.has_tag("zD") for rd in reads.all) == 0

    dp = DupmarkPreprocessor(
        prefilter_params=PrefilterParamsDupmark(min_mapq=0, min_avg_clipq=0, min_baseq=0),
        fixed_params=FixedParamsDupmark(duplication_window_size=boundary_wobble),
    )
    dp._var_params = unsafe_construct_params(RunParamsDupmark, reads=reads)

    dp.process()
    assert sum(rd.has_tag("zD") for rd in reads.all) == 1
