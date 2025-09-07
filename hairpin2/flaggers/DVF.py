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
from pysam import AlignedSegment
from typing import override
from enum import IntEnum, auto
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.abstractflaggers import FixedParams, FlagResult, ReadAwareProcess, read_tagger, require_read_properties, variant_flagger
from hairpin2 import ref2seq as r2s
from typing import cast

STUTTER_DUP_TAG = 'zD'
_FLAG_NAME = "DVF"


class FixedParamsDupmark(FixedParams):
    duplication_window_size: int = 6  # -1 disables


def tag_dups(
    run_params: RunParamsShared,
    fixed_params: FixedParamsDupmark
):
    if not len(run_params.reads.all) > 1:
        return

    for reads in run_params.reads.values():
        if not len(reads) > 1:
            continue

        # prep data
        sample_pair_endpoints: list[tuple[AlignedSegment, tuple[int, int, int, int]]] = []
        for read in run_params.reads.all:
            next_ref_end = r2s.ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
            pair_endpoints: tuple[int, int, int, int] = cast(
            tuple[int, int, int, int],
                tuple(
                    sorted([
                        read.reference_start,
                        cast(int, read.reference_end),
                        read.next_reference_start,
                        next_ref_end
                    ])
                )
            )
            sample_pair_endpoints.append((read, pair_endpoints))
        sample_pair_endpoints = sorted(sample_pair_endpoints, key=lambda x: x[1][0])

        # test data
        # NOTE: ideally keep highest quality read from dups, not currently done
        # NOTE: do you want to only look at supporting reads for this dupmarking? Discuss
        dup_endpoint_comparison_pool: list[tuple[int, int, int, int]] = []
        dup_endpoint_comparison_pool.append(sample_pair_endpoints[0][1])  # endpoints
        for i in range(1, len(sample_pair_endpoints)):
            read = sample_pair_endpoints[i][0]
            testing_endpoints: tuple[int, int, int, int] = sample_pair_endpoints[i][1]
            max_diff_per_comparison: list[int] = []
            for comparison_endpoints in dup_endpoint_comparison_pool:
                endpoint_diffs: tuple[int, int, int, int] = cast(
                    tuple[int, int, int, int],
                    tuple(
                        [abs(x - y) for x, y in zip(comparison_endpoints, testing_endpoints)]
                    )
                )
                assert len(endpoint_diffs) == 4  # sanity check
                max_diff_per_comparison.append(max(endpoint_diffs))
            if all([x <= fixed_params.duplication_window_size for x in max_diff_per_comparison]):
                # then the read pair being examined is a duplicate of the others in the pool
                dup_endpoint_comparison_pool.append(testing_endpoints)  # add it to the comparison pool
                read.set_tag(STUTTER_DUP_TAG, 1, 'i')  # TODO: surface this a bit more - this is our stutter dup tag
                # read.is_duplicate = True  # dupmark
            else:
                # read at i is not dup of reads in dup_endpoint_test_pool
                # start again, test read at i against reads subsequent from i in ends_sorted
                dup_endpoint_comparison_pool = [testing_endpoints]


# TODO: make end user not need to import and inherit from IntEnum, provide some kind of construction method?
class CodesDVF(IntEnum):
    INSUFFICIENT_READS = 0
    DUPLICATION = auto()


class ResultDVF(
    FlagResult,
    flag_name=_FLAG_NAME,
    result_codes=tuple(CodesDVF)
):
    alt: str
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.loss_ratio}"  # TODO: report which samples?


class FixedParamsDVF(FixedParams):
    # NOTE: neither of these options prevent read removal due to duplication,
    # so `DVF.test()` always functions as QC and may still drop reads
    # - I think this is fine, just document more
    read_loss_threshold: float = 0.49  # percent threshold of N duplicate reads compared to N input reads for a given variant and sample, above which we call DVF
    nsamples_threshold: int = 0  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


# detect PCR duplicates previously missed due to slippage
# this implementation assumes that sorting on first element of each sublist
# is appropriate, per Peter's initial implementation.
# is an all against all comparison between all read lists more appropriate?
# and between pairs of readlists, why is comparing sorted pairs most appropriate?
def test_duplicated_support_frac(
    run_params: RunParamsShared,
    fixed_params: FixedParamsDVF
) -> ResultDVF:
    """
    A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.
    """
    # NOTE: surely this shouldn't be across samples... I don't know, maybe?
    nsamples_with_duplication = 0
    nreads_by_sample: dict[str, int] = { k: len(v) for k, v in run_params.reads.items() }
    loss_ratio: list[float] = []

    if not any([nreads > 1 for nreads in nreads_by_sample.values()]):
        fresult = ResultDVF(
            flag=None,
            code=CodesDVF.INSUFFICIENT_READS,
            alt=run_params.alt,
            loss_ratio=0
        )
    else:
        code = CodesDVF.DUPLICATION  # testing possible, and this is the only relevant code
        for sample_key, reads in run_params.reads.items():
            ntotal = nreads_by_sample[sample_key]
            sample_loss_ratio = 0
            if ntotal > 1:
                ndup = sum(read.has_tag('zD') for read in reads)
                ntrue = abs(ndup - ntotal)
                assert ntotal > ndup > -1  # sanity check - TODO: should probably throw an interpretable error
                assert ntotal >= ntrue > -1
                sample_loss_ratio = ndup / ntotal
                if sample_loss_ratio > fixed_params.read_loss_threshold or ntrue < 2:
                    nsamples_with_duplication += 1

            loss_ratio.append(sample_loss_ratio)

        if nsamples_with_duplication > fixed_params.nsamples_threshold:
            flag = True
        else:
            flag = False
        fresult = ResultDVF(
            flag=flag,
            code=code,
            alt=run_params.alt,
            loss_ratio=sum(loss_ratio) / len(loss_ratio)  # TODO: discuss whether averaging is the best choice
        )

    return fresult


# (old note but worth keeping)
# NOTE: since the tagger and flagger are tied to the same class, can't currently specify different/requires excludes
# this is not ideal since we probably do want to dupmark low quality reads - though it does match the previous implementation
# since low qual was already removed before dupmarking.
# The solution seems to be to subsume required_read_properties into the other decorators, so each decorator gives independent
# reqs. Oh, but then you can't solve the run order in the same way... I guess they should be separate.
# Oh well, do this way first to check for regressions. And it may be possible still to solve run order just with a slightly
# different approach in the internals - which would be nice because then you can keep stuff conceptually tied
# TODO/BUG: overall params mapping in ReadAwareProcess doesn't handle namespacing of params with the same name!!
# NOTE: used to exclude low quality - now not doing that per Phuong - NO WAIT YES I AM - have to check no regressions


# TODO/BUG/NOTE: excluding zQ, low qual, because with a bad MC this will fail
# but that's a specific sub property of lq which should itself be surfaced
# it's basically whether the read in question properly paired or not
@require_read_properties(require_tags=['MC', 'zS'], exclude_tags=['zQ'])  # require support  # NOTE: now different to below, could introduce regression
@read_tagger(tagger_param_class=FixedParamsDupmark, read_modifier_func=tag_dups, adds_tag=STUTTER_DUP_TAG)
class TaggerDupmark(
    ReadAwareProcess
): pass


# NOTE/BUG/TODO: excludes low qual tag, assigned by FlaggerLQF in it's read tagger method
# which also has a flagger method that relies on dup tag, assigned above.
# Somewhat messy interdependence
@require_read_properties(require_tags=['MC', 'zS'], exclude_tags=['zO', 'zQ'])  # require support, exclude overlapping second member
@variant_flagger(flag_name=_FLAG_NAME, flagger_param_class=FixedParamsDVF, flagger_func=test_duplicated_support_frac, result_type=ResultDVF)
class FlaggerDVF(
    ReadAwareProcess
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

