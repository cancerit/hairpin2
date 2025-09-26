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
from htsflow.configure_funcs import make_read_processor, make_variant_flagger
from htsflow.process import ReadAwareProcess
from htsflow.process_params import FixedParams
from htsflow.structures import FlagResult
from hairpin2.const import FlaggerNamespaces, Tags, TaggerNamespaces
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.sci_funcs import PropConds, tag_reads_as_duplicates, test_proportion_with_tag



# DUPMARK READ PROCESSOR

class FixedParamsDupmark(FixedParams):
    duplication_window_size: int = 6  # -1 disables


def tag_dups(
    run_params: RunParamsShared,
    fixed_params: FixedParamsDupmark
):
    tag_reads_as_duplicates(
        run_params.reads.all,
        fixed_params.duplication_window_size
    )


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_STUTTER_DUP,
    tagger_param_class=FixedParamsDupmark,
    read_modifier_func=tag_dups,
    adds_marks=[Tags.STUTTER_DUP_TAG],
)
class TaggerDupmark(
    ReadAwareProcess,
): pass





# VARIANT FROM DUPLICATION FLAGGER

@dataclass(frozen=True)
class ResultDVF(
    FlagResult,
    flag_name=FlaggerNamespaces.DUPLICATION,
    info_enum=PropConds
):
    alt: str
    reads_seen: int
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        info_bits = hex(self.info_flag.value if self.info_flag is not None else 0)
        return f"{self.alt}|{self.variant_flagged.value}|{info_bits}|{self.reads_seen}|{round(self.loss_ratio, 3)}"


class FixedParamsDVF(FixedParams):
    """
        read_loss_threshold - percent threshold of N lq reads compared to N input reads for a given variant and sample, above which we flag DVF
        min_pass_reads - the absolute mininum number of reads required for a variant not to be flagged DVF
    """
    read_loss_threshold: float
    min_pass_reads: int
    nsamples_threshold: int  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


# detect PCR duplicates previously missed due to slippage
# this implementation assumes that sorting on first element of each sublist
# is appropriate, per Peter's initial implementation.
# is an all against all comparison between all read lists more appropriate?
# and between pairs of readlists, why is comparing sorted pairs most appropriate?
# def test_duplicated_support_frac(
#     run_params: RunParamsShared,
#     fixed_params: FixedParamsDVF
# ) -> ResultDVF:
#     # NOTE: surely this shouldn't be across samples... I don't know, maybe?

#     if not run_params.reads.all:
#         fresult = ResultDVF(
#             variant_flagged=TO.NA,
#             info_flag=InfoFlagsDVF.INSUFFICIENT_READS,
#             alt=run_params.alt,
#             reads_seen=0,
#             loss_ratio=0
#         )
#     else:
#         info_bits = InfoFlagsDVF.NODATA
#         outcome = TO.VARIANT_PASS
#         # code = InfoFlagsLQF.LOW_QUAL  # testing possible, and this is the only relevant code
#         reads = run_params.reads.all
#         ntotal = len(reads)
#         ndup = sum(read_has_mark(read, Tags.STUTTER_DUP_TAG) for read in reads)
#         ntrue = ntotal - ndup
#         loss_ratio= ndup / ntotal
#         if loss_ratio > fixed_params.read_loss_threshold:
#             info_bits |= InfoFlagsDVF.THRESHOLD
#             outcome = TO.VARIANT_FAIL
#         if ntrue < fixed_params.min_pass_reads:
#             info_bits |= InfoFlagsDVF.MIN_PASS
#             outcome = TO.VARIANT_FAIL

#         if outcome == TO.VARIANT_PASS:
#             info_bits = ~InfoFlagsDVF(0)

#         fresult = ResultDVF(
#             variant_flagged=outcome,
#             info_flag=info_bits,
#             alt=run_params.alt,
#             reads_seen=0,
#             loss_ratio=loss_ratio 
#         )

#     return fresult


def test_DVF(
    run_params: RunParamsShared,
    fixed_params: FixedParamsDVF
):
    result = test_proportion_with_tag(
        run_params.reads.all,
        [Tags.STUTTER_DUP_TAG],
        fixed_params.read_loss_threshold,
        fixed_params.min_pass_reads
    )

    flag = ResultDVF(
        result.outcome,
        result.reason,
        run_params.alt,
        len(run_params.reads.all),
        result.prop_loss
    )

    return flag


@make_variant_flagger(
    process_namespace=FlaggerNamespaces.DUPLICATION,
    flagger_param_class=FixedParamsDVF,
    flagger_func=test_DVF,
    result_type=ResultDVF,
)
class FlaggerDVF(
    ReadAwareProcess,
):
    """
    duplication variant filter - a portion of the reads supporting the variant
    are suspected to arise from duplicated reads that have escaped dupmarking.

    In regions of low complexity, short repeats and homopolymer tracts can cause PCR stuttering.
    Leading to, for example, an additional A on the read when amplifying a tract of As.
    If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference
    between reads that are in fact duplicates. Because of this, these duplicates then evade dupmarking and give rise to
    spurious variants when calling.

    `min_boundary_deviation` sets the minimum deviation start/end coordinates, above which reads are assumed not to be duplicated
    `read_number_difference_threshold` sets the the threshold for absolute difference between the number of reads supporting the variant
    with and without duplicates removed. If this threshold is exceeded, the flag will be set.
    """
    pass

