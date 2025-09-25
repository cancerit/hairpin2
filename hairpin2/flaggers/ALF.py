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

from dataclasses import dataclass
from typing import override
from enum import Flag, auto
from statistics import median
from hairpin2.abstractions.configure_funcs import make_variant_flagger
from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.process_params import FixedParams
from hairpin2.abstractions.structures import FlagResult
from hairpin2.const import FlaggerNamespaces, TestOutcomes as TO
from hairpin2.flaggers.shared import RunParamsShared

# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


class InfoFlagsALF(Flag):
    INSUFFICIENT_READS = 1
    INSUFFICIENT_AS_TAGS = auto()
    ON_THRESHOLD = auto()


@dataclass(frozen=True)
class ResultALF(
    FlagResult,
    flag_name="ALF",
    info_enum=InfoFlagsALF
):
    alt: str
    avg_as: float | None

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.variant_flagged}|{self.info_flag}|{self.avg_as}"


class FixedParamsALF(FixedParams):
    avg_AS_threshold: float

def test_alignment_score(  # test supporting reads
    run_params: RunParamsShared,
    fixed_params: FixedParamsALF
):
    if len(run_params.reads.all) < 1:
        code = InfoFlagsALF.INSUFFICIENT_READS
        fresult = ResultALF(
            variant_flagged=TO.NA,
            info_flag=code,
            alt=run_params.alt,
            avg_as=None
        )
    else:
        aln_scores: list[float] = []

        for read in run_params.reads.all:
            try:
                aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
            except KeyError:
                pass
        if len(aln_scores) != 0:
            avg_as = median(aln_scores)
            code = InfoFlagsALF.ON_THRESHOLD
            flag = TO.VARIANT_PASS
            if avg_as <= fixed_params.avg_AS_threshold:
                flag = TO.VARIANT_FAIL
            fresult = ResultALF(
                variant_flagged=flag,
                info_flag=code,
                alt=run_params.alt,
                avg_as=avg_as
            )
        else:
            code = InfoFlagsALF.INSUFFICIENT_AS_TAGS
            fresult = ResultALF(
                variant_flagged=TO.NA,
                info_flag=code,
                alt=run_params.alt,
                avg_as=None
            )

    return fresult


@make_variant_flagger(
    process_namespace=FlaggerNamespaces.POOR_ALIGNMENT_SCORE,
    flagger_param_class=FixedParamsALF,
    flagger_func=test_alignment_score,result_type=ResultALF
)
class FlaggerALF(
    ReadAwareProcess,
):
    """
    Alignment score filter based on AS tag
    """
