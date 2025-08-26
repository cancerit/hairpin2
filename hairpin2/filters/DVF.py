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
import hairpin2.abstractfilters as haf
from hairpin2 import ref2seq as r2s
from pydantic.dataclasses import dataclass
from typing import ClassVar, override, cast
from pysam import AlignedSegment
from enum import IntEnum, auto


class DVCodes(IntEnum):
    INSUFFICIENT_READS = 0
    DUPLICATION = auto()


class Result(haf.FilterResult[DVCodes]):
    Name: ClassVar[str] = 'DVF'
    alt: str
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.loss_ratio}"  # TODO: report which samples?


@dataclass(slots=True, frozen=True)
class Params(haf.FilterParams):
    duplication_window_size: int = 6  # -1 disables
    # NOTE: neither of these options prevent read removal due to duplication,
    # so `DVF.test()` always functions as QC and may still drop reads
    # - I think this is fine, just document more
    read_loss_threshold: float = 0.49  # percent threshold of N duplicate reads compared to N input reads for a given variant and sample, above which we call DVF
    nsamples_threshold: int = 0  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


class Filter(haf.FilterTester[dict[str, list[AlignedSegment]], Params, Result]):
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
    # detect PCR duplicates previously missed due to slippage
    # this implementation assumes that sorting on first element of each sublist
    # is appropriate, per Peter's initial implementation.
    # is an all against all comparison between all read lists more appropriate?
    # and between pairs of readlists, why is comparing sorted pairs most appropriate?
    @override
    def test(
        self,
        alt: str,
        variant_reads_by_sample: dict[str, list[AlignedSegment]]
    ) -> tuple[dict[str, list[AlignedSegment]], Result]:
        """
        A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.
        """
        nsamples_with_duplication = 0
        nreads_by_sample: dict[str, int] = { k: len(v) for k, v in variant_reads_by_sample.items() }
        sanitised_reads_by_sample: dict[str, list[AlignedSegment]] = {}
        loss_ratio: list[float] = []

        if not any([nreads > 1 for nreads in nreads_by_sample.values()]):
            fresult = Result(
                flag=None,
                code=DVCodes.INSUFFICIENT_READS,
                alt=alt,
                loss_ratio=0
            )
            sanitised_reads_by_sample = variant_reads_by_sample
        else:
            code = DVCodes.DUPLICATION  # testing possible, and this is the only relevant code
            for sample_key, reads in variant_reads_by_sample.items():
                sample_loss_ratio = 0
                if nreads_by_sample[sample_key] > 1:
                    # prep data
                    sample_pair_endpoints: list[tuple[int, tuple[int, int, int, int]]] = []
                    for idx, read in enumerate(reads):
                        next_ref_end = r2s.ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
                        if read.reference_end is None:
                            raise ValueError(f"read {read.query_name} has no reference end coordinate. all reads must be fully described.")
                        else:
                            pair_endpoints: tuple[int, int, int, int] = cast(
                            tuple[int, int, int, int],
                                tuple(
                                    sorted([
                                        read.reference_start,
                                        read.reference_end,
                                        read.next_reference_start,
                                        next_ref_end
                                    ])
                                )
                            )
                            sample_pair_endpoints.append((idx, pair_endpoints))
                    sample_pair_endpoints = sorted(sample_pair_endpoints, key=lambda x: x[1][0])

                    # test data
                    dup_idcs: list[int] = []
                    dup_endpoint_test_pool: list[tuple[int, int, int, int]] = []
                    dup_endpoint_test_pool.append(sample_pair_endpoints[0][1])
                    for i in range(1, len(sample_pair_endpoints)):
                        original_index = sample_pair_endpoints[i][0]
                        testing_endpoints: tuple[int, int, int, int] = sample_pair_endpoints[i][1]
                        max_diff_per_comparison: list[int] = []
                        for comparison_endpoints in dup_endpoint_test_pool:
                            endpoint_diffs: tuple[int, int, int, int] = cast(
                                tuple[int, int, int, int],
                                tuple(
                                    [abs(x - y) for x, y in zip(comparison_endpoints, testing_endpoints)]
                                )
                            )
                            max_diff_per_comparison.append(max(endpoint_diffs))
                        if all([x <= self.fixed_params.duplication_window_size for x in max_diff_per_comparison]):
                            # then the read pair being examined is a duplicate of the others in the pool
                            dup_endpoint_test_pool.append(testing_endpoints)
                            dup_idcs.append(original_index)  # store original index of this duplicate
                        else:
                            # read at i is not dup of reads in dup_endpoint_test_pool
                            # start again, test read at i against reads subsequent from i in ends_sorted
                            dup_endpoint_test_pool = [testing_endpoints]
                    sanitised_reads = [read for i, read in enumerate(reads) if i not in dup_idcs]
                    n_loss = nreads_by_sample[sample_key] - len(sanitised_reads)
                    assert nreads_by_sample[sample_key] > n_loss > -1
                    sample_loss_ratio = n_loss / nreads_by_sample[sample_key]
                    if sample_loss_ratio > self.fixed_params.read_loss_threshold or len(sanitised_reads) == 1:
                        nsamples_with_duplication += 1
                else:
                    sanitised_reads = reads

                sanitised_reads_by_sample[sample_key] = sanitised_reads
                loss_ratio.append(sample_loss_ratio)

            if nsamples_with_duplication > self.fixed_params.nsamples_threshold:
                flag = True
            else:
                flag = False
            fresult = Result(
                flag=flag,
                code=code,
                alt=alt,
                loss_ratio=sum(loss_ratio) / len(loss_ratio)  # TODO: discuss whether averaging is the best choice
            )

        return sanitised_reads_by_sample, fresult
