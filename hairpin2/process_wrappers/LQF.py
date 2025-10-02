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
from typing import Annotated, override

from pydantic import BeforeValidator

from hairpin2.const import FlaggerNamespaces, TaggerNamespaces, Tags
from hairpin2.infrastructure.configure_funcs import make_read_processor, make_variant_flagger
from hairpin2.infrastructure.constraints import bound
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.infrastructure.process_params import FixedParams
from hairpin2.infrastructure.structures import FlagResult
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import ProportionBasedTest, TagLowQualReads

# LOW QUAL READ TAGGER


class QualParams(FixedParams):
    min_mapping_quality: Annotated[int, BeforeValidator(lambda x: bound(x, 0, 255))]
    min_avg_clip_quality: Annotated[int, BeforeValidator(lambda x: bound(x, 0, 40))]
    min_base_quality: Annotated[int, BeforeValidator(lambda x: bound(x, 0, 40))]


def tag_lq(
    run_params: RunParamsShared,
    params: QualParams,  # placeholder maybe
):
    for read in run_params.reads.all:
        _ = TagLowQualReads.check_low_qual_read(
            read,
            run_params.record.start,
            run_params.alt,
            run_params.mut_type,
            params.min_base_quality,
            params.min_mapping_quality,
            params.min_avg_clip_quality,
            mark=True,
        )


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_LOW_QUAL,
    tagger_param_class=QualParams,
    read_modifier_func=tag_lq,
    adds_marks=[Tags.LOW_QUAL_TAG],
)
class TaggerLowQual(ReadAwareProcess):
    pass


# LOW QUAL VARIANT FLAGGER


@dataclass(frozen=True)
class ResultLQF(
    FlagResult, flag_name=FlaggerNamespaces.LOW_QUAL, info_enum=ProportionBasedTest.ResultPack.Info
):
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self, alt: str) -> str:
        info_bits = hex(self.info_flag.value if self.info_flag is not None else 0)
        return f"{alt}|{self.variant_flagged.value}|{info_bits}|{self.reads_seen}|{round(self.loss_ratio, 3)}"


class FixedParamsLQF(FixedParams):
    """
    read_loss_threshold - percent threshold of N lq reads compared to N input reads for a given variant and sample, above which we flag LQF
    min_pass_reads - the absolute mininum number of reads required for a variant not to be flagged LQF
    """

    read_loss_threshold: Annotated[float, BeforeValidator(lambda x: bound(x, 0.0, 1.0))]
    min_pass_reads: Annotated[int, BeforeValidator(lambda x: bound(x, 0))]
    nsamples_threshold: int  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


def test_variant_LQF(run_params: RunParamsShared, fixed_params: FixedParamsLQF):
    result = ProportionBasedTest.test_variant_reads(
        run_params.reads.all,
        [Tags.LOW_QUAL_TAG, Tags.STUTTER_DUP_TAG],
        fixed_params.read_loss_threshold,
        fixed_params.min_pass_reads,
    )

    flag = ResultLQF(
        result.outcome, result.reason, len(run_params.reads.all), result.prop_loss
    )

    return flag


# exclude overlapping second pair member, require support  NOTE: do you actually want to exclude overlap when assessing this? it's not quite double counting per se if the overlapping read does show support...
@make_variant_flagger(
    process_namespace=FlaggerNamespaces.LOW_QUAL,
    flagger_func=test_variant_LQF,
    flagger_param_class=FixedParamsLQF,
    result_type=ResultLQF,
)
class FlaggerLQF(ReadAwareProcess):
    pass
