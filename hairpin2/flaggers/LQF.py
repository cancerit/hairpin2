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
from statistics import mean
from pysam import AlignedSegment
import hairpin2.abstractflaggers as haf
from typing import Any, cast, override
from enum import IntEnum, auto
from hairpin2.flaggers.shared import PrefilterParamsShared, RunParamsShared
from hairpin2.readqc import ValidatorFlags
from hairpin2.ref2seq import ref2querypos, ref_end_via_cigar


_FLAG_NAME = 'LQF'


# TODO: make end user not need to import and inherit from IntEnum, provide some kind of construction method?
class CodesLQF(IntEnum):
    INSUFFICIENT_READS = 0
    LOW_QUAL = auto()


class ResultLQF(
    haf.FlagResult,
    flag_name=_FLAG_NAME,
    result_codes=tuple(CodesLQF)
):
    alt: str
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.loss_ratio}"  # TODO: report which samples?


class FixedParamsLQF(haf.FixedParams):
    read_loss_threshold: float = 0.49  # percent threshold of N lq reads compared to N input reads for a given variant and sample, above which we call DVF
    nsamples_threshold: int = 0  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


def qc_read(
    read: AlignedSegment,
    vcf_start: int,
    alt: str,
    mut_type: str,
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

    """
    When testing a variant, test a read for various features specific to the variant
    at hand which would identify that read as a poor source of support for that
    variant. Used to disqualify reads for use in testing a variant.
    """
    try:
        mate_cig = str(read.get_tag('MC'))
    except KeyError:
        mate_cig = None
    if any(flg is None for flg in
            [read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig]):
        invalid_flag |= ValidatorFlags.READ_FIELDS_MISSING
    else:
        if not (read.flag & 0x2) or read.flag & 0xE00:
            invalid_flag |= ValidatorFlags.FLAG

        if read.mapping_quality < min_mapqual:
            invalid_flag |= ValidatorFlags.MAPQUAL

        if ('S' in read.cigarstring and  # pyright: ignore[reportOperatorIssue]
                mean(read.query_alignment_qualities) < min_clipqual):  # pyright: ignore[reportUnknownMemberType, reportArgumentType]
            invalid_flag |= ValidatorFlags.CLIPQUAL

        # avoid analysing both read1 and mate if they both cover the variant
        # NOTE: introduces strand bias!!
        if (not (invalid_flag & ValidatorFlags.FLAG)
                and not (read.flag & 0x40)):
            read_range = range(read.reference_start,
                               read.reference_end)  # pyright: ignore[reportArgumentType]
            mate_range = range(read.next_reference_start,
                               ref_end_via_cigar(mate_cig,  # pyright: ignore[reportArgumentType]
                                                     read.next_reference_start))
            ref_overlap = set(read_range).intersection(mate_range)
            if vcf_start in ref_overlap:
                invalid_flag |= ValidatorFlags.OVERLAP


    if mut_type == 'S':
        try:
            mut_pos = ref2querypos(read, vcf_start)
        except ValueError:
            invalid_flag |= ValidatorFlags.NOT_ALIGNED
        else:
            if any(
                [bq < min_basequal
                for bq
                in cast(array[Any], read.query_qualities)[mut_pos:mut_pos + len(alt)]]
            ):
                invalid_flag |= ValidatorFlags.BASEQUAL

    return invalid_flag


LOW_QUAL_TAG = 'zQ'


def tag_lq(
    run_params: RunParamsShared,
    params: PrefilterParamsShared # placeholder maybe
):
    for read in run_params.reads.all:
        if qc_read(
            read,
            run_params.record.start,
            run_params.alt,
            run_params.mut_type,
            params.min_base_quality,
            params.min_mapping_quality,
            params.min_avg_clip_quality,
            
        ) != ValidatorFlags.CLEAR:  # if bad
            read.set_tag(LOW_QUAL_TAG, 1, 'i')



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
            flag=None,
            code=CodesLQF.INSUFFICIENT_READS,
            alt=run_params.alt,
            loss_ratio=0
        )
    else:
        code = CodesLQF.LOW_QUAL  # testing possible, and this is the only relevant code
        for sample_key, reads in run_params.reads.items():
            ntotal = nreads_by_sample[sample_key]
            sample_loss_ratio = 0
            if ntotal > 1:
                nlq = sum((read.has_tag(LOW_QUAL_TAG) or read.has_tag('zD')) for read in reads)
                ntrue = abs(nlq - ntotal)
                assert ntotal >= nlq > -1  # sanity check - TODO: should probably throw an interpretable error
                assert ntotal >= ntrue > -1
                assert ntotal == nlq + ntrue
                sample_loss_ratio = nlq / ntotal
                if sample_loss_ratio > fixed_params.read_loss_threshold or ntrue < 2:
                    nsamples_with_lowqual += 1

            loss_ratio.append(sample_loss_ratio)

        if nsamples_with_lowqual > fixed_params.nsamples_threshold:
            flag = True
        else:
            flag = False
        fresult = ResultLQF(
            flag=flag,
            code=code,
            alt=run_params.alt,
            loss_ratio=sum(loss_ratio) / len(loss_ratio)  # TODO: discuss whether averaging is the best choice
        )

    return fresult

@haf.require_read_properties(require_tags=['zS'])  # require support
@haf.read_tagger(tagger_param_class=PrefilterParamsShared, read_modifier_func=tag_lq, adds_tag=LOW_QUAL_TAG)  # placeholder
class TaggerLowQual(
    haf.ReadAwareProcess,
    process_name="mark-low-qual"
): pass


@haf.require_read_properties(require_tags=['zS'], exclude_tags=['zO'])  # exclude overlapping second pair member, require support  NOTE: do you actually want to exclude overlap when assessing this? it's not quite double counting per se if the overlapping read does show support...
@haf.variant_flagger(
    flag_name=_FLAG_NAME,
    flagger_func=test_variant_LQF,
    flagger_param_class=FixedParamsLQF,
    result_type=ResultLQF
)
class FlaggerLQF(
    haf.ReadAwareProcess,
    process_name=_FLAG_NAME
):
    """
    """
    # NOTE: checks dups, meaning dupmarking must have run first - not at all surfaced in this implementation
