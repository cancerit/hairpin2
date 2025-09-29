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

from hairpin2.const import FlaggerNamespaces
from hairpin2.infrastructure.configure_funcs import make_variant_flagger
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.infrastructure.process_params import FixedParams
from hairpin2.infrastructure.structures import FlagResult
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import AlignmentScoreTest


@dataclass(frozen=True)
class ResultALF(
    FlagResult,
    flag_name="ALF",
    info_enum=AlignmentScoreTest.ResultPack.Info,
):
    alt: str
    reads_seen: int
    avg_as: float | None

    @override
    def getinfo(self) -> str:
        info_bits = hex(self.info_flag.value if self.info_flag is not None else 0)
        avg_as = round(self.avg_as, 3) if self.avg_as else "NA"
        return f"{self.alt}|{self.variant_flagged}|{info_bits}|{self.reads_seen}|{avg_as}"


class FixedParamsALF(FixedParams):
    avg_AS_threshold: float


def test_ALF(  # test supporting reads
    run_params: RunParamsShared, fixed_params: FixedParamsALF
):
    result = AlignmentScoreTest.test_variant_reads(
        run_params.reads.all, fixed_params.avg_AS_threshold
    )

    flag = ResultALF(
        result.outcome, result.reason, run_params.alt, len(run_params.reads), result.avg_as
    )

    return flag


@make_variant_flagger(
    process_namespace=FlaggerNamespaces.POOR_ALIGNMENT_SCORE,
    flagger_param_class=FixedParamsALF,
    flagger_func=test_ALF,
    result_type=ResultALF,
)
class FlaggerALF(
    ReadAwareProcess,
):
    """
    Alignment score filter based on AS tag
    """
