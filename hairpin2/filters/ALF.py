# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
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
import hairpin2.abstractfilters as haf
from pysam import AlignedSegment
from pydantic.dataclasses import dataclass
from typing import override, ClassVar
from collections.abc import Sequence
from enum import IntEnum, auto
from statistics import median
# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


class ALCodes(IntEnum):
    INSUFFICIENT_READS = 0
    INSUFFICIENT_AS_TAGS = auto()
    ON_THRESHOLD = auto()


class Result(haf.FilterResult[ALCodes]):
    Name: ClassVar[str] = 'ALF'
    alt: str
    avg_as: float | None

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.avg_as}"


@dataclass(frozen=True, slots=True)
class Params(haf.FilterParams):
    al_thresh: float = 0.93


class Filter(haf.FilterTester[Sequence[AlignedSegment], Params, Result]):
    """
    Describe filter
    """
    @override
    def test[T: Sequence[AlignedSegment]](
        self,
        alt: str,
        variant_reads: T,
    ) -> tuple[T, Result]:
        if len(variant_reads) < 1:
            code = ALCodes.INSUFFICIENT_READS
            fresult = Result(
                flag=None,
                code=code,
                alt=alt,
                avg_as=None
            )
        else:
            aln_scores: list[float] = []

            for read in variant_reads:
                try:
                    aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
                except KeyError:
                    pass
            if len(aln_scores) != 0:
                avg_as = median(aln_scores)
                code = ALCodes.ON_THRESHOLD
                flag = False
                if avg_as <= self.fixed_params.al_thresh:
                    flag = True
                fresult = Result(
                    flag=flag,
                    code=code,
                    alt=alt,
                    avg_as=avg_as
                )
            else:
                code = ALCodes.INSUFFICIENT_AS_TAGS
                flag = None
                fresult = Result(
                    flag=flag,
                    code=code,
                    alt=alt,
                    avg_as=None
                )

        return variant_reads, fresult
