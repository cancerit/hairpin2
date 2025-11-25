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
import copy

from hairpin2.infrastructure.structures import ExtendedRead, TestOutcomes
from hairpin2.sci_funcs import AlignmentScoreTest


class TestAlignmentScoreTest:
    def test_low_alignment_score(self, read: ExtendedRead):
        test_read = copy.deepcopy(read)
        test_read.set_tag("AS", 50)  # low AS
        reads_in = [test_read, test_read]
        result = AlignmentScoreTest.test_variant_reads(reads_in, 0.93)
        assert result.outcome == TestOutcomes.VARIANT_FAIL
        assert result.reason & AlignmentScoreTest.ResultPack.Info.ON_THRESHOLD


    def test_high_alignment_score(self, read: ExtendedRead):
        test_read = copy.deepcopy(read)
        test_read.set_tag("AS", 99)  # low AS
        reads_in = [test_read, test_read]
        result = AlignmentScoreTest.test_variant_reads(reads_in, 0.93)
        assert result.outcome == TestOutcomes.VARIANT_PASS
        assert result.reason & AlignmentScoreTest.ResultPack.Info.ON_THRESHOLD


    def test_insufficient_alignment_score(self, read: ExtendedRead):
        test_read = copy.deepcopy(read)
        reads_in = [test_read, test_read]
        result = AlignmentScoreTest.test_variant_reads(reads_in, 0.93)
        assert result.outcome == TestOutcomes.NA
        assert result.reason & AlignmentScoreTest.ResultPack.Info.INSUFFICIENT_AS_TAGS


    def test_insufficient_reads(self):
        reads_in = []
        result = AlignmentScoreTest.test_variant_reads(reads_in, 0.93)
        assert result.outcome == TestOutcomes.NA
        assert result.reason & AlignmentScoreTest.ResultPack.Info.INSUFFICIENT_READS
