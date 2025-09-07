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
# TODO: minimise necessary imports for end user...
from typing import override
from enum import IntEnum, auto
from statistics import median
from hairpin2.flaggers.shared import RunParamsShared

# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


_FLAG_NAME = "ALF"


class CodesALF(IntEnum):
    INSUFFICIENT_READS = 0
    INSUFFICIENT_AS_TAGS = auto()
    ON_THRESHOLD = auto()


class ResultALF(
    haf.FlagResult,
    flag_name="ALF",
    result_codes=tuple(CodesALF)
):
    alt: str
    avg_as: float | None

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.avg_as}"


class FixedParamsALF(haf.FixedParams):
    al_thresh: float = 0.93

def test_alignment_score(  # test supporting reads
    run_params: RunParamsShared,
    fixed_params: FixedParamsALF
):
    if len(run_params.reads.all) < 1:
        code = CodesALF.INSUFFICIENT_READS
        fresult = ResultALF(
            flag=None,
            code=code,
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
            code = CodesALF.ON_THRESHOLD
            flag = False
            if avg_as <= fixed_params.al_thresh:
                flag = True
            fresult = ResultALF(
                flag=flag,
                code=code,
                alt=run_params.alt,
                avg_as=avg_as
            )
        else:
            code = CodesALF.INSUFFICIENT_AS_TAGS
            flag = None
            fresult = ResultALF(
                flag=flag,
                code=code,
                alt=run_params.alt,
                avg_as=None
            )

    return fresult


# TODO: make actual mixins not importable so I stop accidentally importing them
# NOTE: consider mapping more informative tags, such as 'SUPPORT', to 2 char htslib tags internally (e.g. 'zS')
@haf.require_read_properties(require_tags=['zS'], exclude_tags=['zD', 'zQ', 'zO'])  # require support, exclude stutter dups, low qual, second overlapping fragment member
@haf.variant_flagger(flag_name=_FLAG_NAME, flagger_param_class=FixedParamsALF, flagger_func=test_alignment_score,result_type=ResultALF)
class FlaggerALF(
    haf.ReadAwareProcess
):
    # TODO: docstring
    """
    Alignment score filter based on AS tag
    """
