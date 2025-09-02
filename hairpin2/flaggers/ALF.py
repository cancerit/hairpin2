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
import hairpin2.abstractflaggers as haf
from pydantic.dataclasses import dataclass
from typing import override, ClassVar
from enum import IntEnum, auto
from statistics import median
from hairpin2.flaggers.shared import AltVarParams

# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


class CodesALF(IntEnum):
    INSUFFICIENT_READS = 0
    INSUFFICIENT_AS_TAGS = auto()
    ON_THRESHOLD = auto()


class ResultALF(haf.FilterResult[CodesALF]):
    Name: ClassVar[str] = 'ALF'
    alt: str
    avg_as: float | None

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.avg_as}"


@dataclass(frozen=True, slots=True)
class FixedParamsALF(haf.FixedParams):
    al_thresh: float = 0.93


class VarParamsALF(AltVarParams): pass


class FlaggerALF(haf.Flagger[FixedParamsALF, VarParamsALF, ResultALF]):
    # TODO: docstring
    """
    Alignment score filter based on AS tag
    """
    @override
    def test(
        self
    ):
        if len(self.var_params.reads.all) < 1:
            code = CodesALF.INSUFFICIENT_READS
            fresult = ResultALF(
                flag=None,
                code=code,
                alt=self.var_params.alt,
                avg_as=None
            )
        else:
            aln_scores: list[float] = []

            for read in self.var_params.reads.all:
                try:
                    aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
                except KeyError:
                    pass
            if len(aln_scores) != 0:
                avg_as = median(aln_scores)
                code = CodesALF.ON_THRESHOLD
                flag = False
                if avg_as <= self.fixed_params.al_thresh:
                    flag = True
                fresult = ResultALF(
                    flag=flag,
                    code=code,
                    alt=self.var_params.alt,
                    avg_as=avg_as
                )
            else:
                code = CodesALF.INSUFFICIENT_AS_TAGS
                flag = None
                fresult = ResultALF(
                    flag=flag,
                    code=code,
                    alt=self.var_params.alt,
                    avg_as=None
                )

        return fresult
