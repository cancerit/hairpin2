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
from array import array
from dataclasses import dataclass
from statistics import mean
from pysam import AlignedSegment
from typing import Any, cast, override
from enum import Flag, auto
from hairpin2.abstractions.configure_funcs import make_read_processor, make_variant_flagger
from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.process_params import FixedParams
from hairpin2.abstractions.structures import ExtendedRead, FlagResult, mark_read, read_has_mark, record_operation
from hairpin2.const import FlaggerNamespaces, MutTypes, Tags, TaggerNamespaces, ValidatorFlags, TestOutcomes as TO
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.utils.ref2seq import ref2querypos


class QualParams(FixedParams):
    min_mapping_quality: int
    min_avg_clip_quality: int
    min_base_quality: int


# TODO: make end user not need to import and inherit from IntEnum, provide some kind of construction method?
class InfoFlagsLQF(Flag):
    INSUFFICIENT_READS = 1
    LOW_QUAL = auto()  # TODO: stop using auto


@dataclass(frozen=True)
class ResultLQF(
    FlagResult,
    flag_name=FlaggerNamespaces.LOW_QUAL,
    info_enum=InfoFlagsLQF
):
    alt: str
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.variant_flagged}|{self.info_flag}|{self.loss_ratio}"  # TODO: report which samples?


class FixedParamsLQF(FixedParams):
    """
        read_loss_threshold - percent threshold of N lq reads compared to N input reads for a given variant and sample, above which we flag LQF
        min_pass_reads - the absolute mininum number of reads required for a variant not to be flagged LQF
    """
    read_loss_threshold: float
    min_pass_reads: int
    nsamples_threshold: int  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


def qc_read(
    read: ExtendedRead | AlignedSegment,
    vcf_start: int,
    alt: str,
    mut_type: MutTypes,
    min_basequal: int,
    min_mapqual: int,
    min_clipqual: int,
) -> ValidatorFlags:
    """
    When testing variants, test a read for various general - i.e. not specific to a
    particular variant - features identifying the read as a poor source of support
    for a variant. Used to disqualify reads for use in testing a variant.
    
    does not check for the following, as pysam already guards against:
        - quality and seq length mismatch
        - reference id is none
    """
    invalid_flag = ValidatorFlags.CLEAR  # 0 - evaluates false


    if not (read.flag & 0x2) or read.flag & 0xE00:
        invalid_flag |= ValidatorFlags.FLAG

    if read.mapping_quality < min_mapqual:
        invalid_flag |= ValidatorFlags.MAPQUAL

    if ('S' in read.cigarstring and  # pyright: ignore[reportOperatorIssue]
            mean(read.query_alignment_qualities) < min_clipqual):  # pyright: ignore[reportUnknownMemberType, reportArgumentType]
        invalid_flag |= ValidatorFlags.CLIPQUAL

    if mut_type == MutTypes.SUB:
        try:
            mut_pos = ref2querypos(read, vcf_start)
        except ValueError:
            raise ValueError(f"read: {read.query_name} appears not to cover variant at {vcf_start} - this function should only be run on reads confirmed to be supporting the variant")
        else:
            if any(
                [bq < min_basequal
                for bq
                in cast(array[Any], read.query_qualities)[mut_pos:mut_pos + len(alt)]]
            ):
                invalid_flag |= ValidatorFlags.BASEQUAL

    return invalid_flag


# TODO: predefined functions for handing back bools from scientist funcs, and I'll handle the tagging
def tag_lq(
    run_params: RunParamsShared,
    params: QualParams # placeholder maybe
):
    for read in run_params.reads.all:
        flag = qc_read(
            read,
            run_params.record.start,
            run_params.alt,
            run_params.mut_type,
            params.min_base_quality,
            params.min_mapping_quality,
            params.min_avg_clip_quality,
            
        )
        if flag != ValidatorFlags.CLEAR:
            mark_read(read, Tags.LOW_QUAL_TAG)
        record_operation(read, 'mark-low-qual')



def test_variant_LQF(
    run_params: RunParamsShared,
    fixed_params: FixedParamsLQF
) -> ResultLQF:
    """
    A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.
    """
    # NOTE: surely this shouldn't be across samples... I don't know, maybe?
    nsamples_with_lowqual = 0
    nreads_by_sample: dict[str, int] = { k: len(v) for k, v in run_params.reads.items() }
    loss_ratio: list[float] = []

    if not any([nreads > 1 for nreads in nreads_by_sample.values()]):
        fresult = ResultLQF(
            variant_flagged=TO.NA,
            info_flag=InfoFlagsLQF.INSUFFICIENT_READS,
            alt=run_params.alt,
            loss_ratio=0
        )
    else:
        code = InfoFlagsLQF.LOW_QUAL  # testing possible, and this is the only relevant code
        for sample_key, reads in run_params.reads.items():
            ntotal = nreads_by_sample[sample_key]
            sample_loss_ratio = 0
            if ntotal > 1:
                nlq = sum((read_has_mark(read, Tags.LOW_QUAL_TAG) or read_has_mark(read, Tags.STUTTER_DUP_TAG)) for read in reads)
                ntrue = abs(nlq - ntotal)
                sample_loss_ratio = nlq / ntotal
                if sample_loss_ratio > fixed_params.read_loss_threshold or ntrue < fixed_params.min_pass_reads:
                    nsamples_with_lowqual += 1

            loss_ratio.append(sample_loss_ratio)

        if nsamples_with_lowqual > fixed_params.nsamples_threshold:
            flag = TO.VARIANT_FAIL
        else:
            flag = TO.VARIANT_PASS
        fresult = ResultLQF(
            variant_flagged=flag,
            info_flag=code,
            alt=run_params.alt,
            loss_ratio=sum(loss_ratio) / len(loss_ratio)  # TODO: discuss whether averaging is the best choice
        )

    return fresult


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_LOW_QUAL,
    tagger_param_class=QualParams,
    read_modifier_func=tag_lq,
    adds_marks=[Tags.LOW_QUAL_TAG],
)
class TaggerLowQual(
    ReadAwareProcess
): pass


# exclude overlapping second pair member, require support  NOTE: do you actually want to exclude overlap when assessing this? it's not quite double counting per se if the overlapping read does show support...
@make_variant_flagger(
    process_namespace=FlaggerNamespaces.LOW_QUAL,
    flagger_func=test_variant_LQF,
    flagger_param_class=FixedParamsLQF,
    result_type=ResultLQF,
)
class FlaggerLQF(
    ReadAwareProcess
):
    """
    """
    # NOTE: checks dups, meaning dupmarking must have run first - not at all surfaced in this implementation
    # hence move to boolean tags, and away from exclude tags
    # DOCSTRINGS
