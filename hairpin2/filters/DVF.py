import hairpin2.abstractfilters as haf
from dataclasses import field
from typing import ClassVar, override, cast
from pysam import AlignedSegment
from hairpin2 import ref2seq as r2s
from enum import IntEnum, auto
# pyright: reportExplicitAny=false
# pyright: reportAny=false
# # pyright: reportUnnecessaryIsInstance=false


class DVCodes(IntEnum):
    INSUFFICIENT_READS = 0
    DUPLICATION = auto()

class DVResult(haf.FilterResult[DVCodes]):
    Codes: ClassVar[type[DVCodes]] = DVCodes
    alt: str
    name: str = field(default='DVF', init=False)

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}"  # TODO: report which samples?

class DVParams(haf.FilterParams):
    min_boundary_deviation: int = 6  # TODO: can be used to turn off - document
    # TODO: n.b. neither of these options prevent read removal due to duplication, so test still functions as QC (and I think this is fine, just document more)
    read_number_difference_threshhold: int = 0 # change in reads! TODO: express as fraction? discuss with Peter/Phuong
    nsamples_threshold: int = 1  # TODO: document

class DVFilter(haf.FilterTester[dict[str, list[AlignedSegment]], DVParams, DVResult]):
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
    ) -> tuple[dict[str, list[AlignedSegment]], DVResult]:
        """
        A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.
        """
        nsamples_with_duplication = 0
        sanitised_reads_by_sample: dict[str, list[AlignedSegment]] = {}

        if not any([len(reads) > 1 for reads in variant_reads_by_sample.values()]):
            fresult = DVResult(
                flag=None,
                code=DVResult.Codes.INSUFFICIENT_READS,
                alt=alt
            )
            sanitised_reads_by_sample = variant_reads_by_sample
        else:
            code = DVResult.Codes.DUPLICATION  # testing possible, and this is the only relevant code
            for sample_key, reads in variant_reads_by_sample.items():
                if len(reads) > 1:
                    # prep data
                    sample_pair_endpoints: list[tuple[int, tuple[int, int, int, int]]] = []
                    for idx, read in enumerate(reads):
                        next_ref_end = r2s.ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
                        if read.reference_end is None:
                            raise ValueError(f"read {read.query_name} has no reference end coordinate. all reads must be fully described.")
                        else:
                            pair_endpoints: tuple[int, int, int, int] = cast(
                                tuple[int, int, int, int],
                                (sorted([
                                    read.reference_start,
                                    read.reference_end,
                                    read.next_reference_start,
                                    next_ref_end
                                ]),)
                            )
                            sample_pair_endpoints.append((idx, pair_endpoints))
                    sample_pair_endpoints = sorted(sample_pair_endpoints, key=lambda x: x[1][0])

                    # test data
                    # TODO: I don't agree with sorting the ends. They should be matched on their field (start, start) - discuss with Phuong/Peter
                    dup_idcs: list[int] = []
                    dup_endpoint_test_pool: list[tuple[int, int, int, int]] = []
                    dup_endpoint_test_pool.append(sample_pair_endpoints[0][1])
                    for i in range(1, len(sample_pair_endpoints)):
                        original_index = sample_pair_endpoints[i][0]
                        testing_endpoints: tuple[int, int, int, int] = sample_pair_endpoints[i][1]
                        max_diff_per_comparison: list[int] = []
                        for comparison_endpoints in dup_endpoint_test_pool:
                            endpoint_diffs: tuple[int, int, int, int] = cast(tuple[int, int, int, int], tuple([abs(x - y) for x, y in zip(comparison_endpoints, testing_endpoints)]))
                            max_diff_per_comparison.append(max(endpoint_diffs))
                        if all([x <= self.fixed_params.min_boundary_deviation for x in max_diff_per_comparison]):
                            # then the read pair being examined is a duplicate of the others in the pool
                            dup_endpoint_test_pool.append(testing_endpoints)
                            dup_idcs.append(original_index)  # store original index of this duplicate
                        else:
                            # read at i is not dup of reads in dup_endpoint_test_pool
                            # start again, test read at i against reads subsequent from i in ends_sorted
                            dup_endpoint_test_pool = [testing_endpoints]
                    sanitised_reads = [read for i, read in enumerate(reads) if i not in dup_idcs]
                    if (len(sanitised_reads) - len(reads)) > self.fixed_params.read_number_difference_threshhold:
                        nsamples_with_duplication += 1
                    sanitised_reads_by_sample[sample_key] = sanitised_reads

            if nsamples_with_duplication > self.fixed_params.nsamples_threshold:
                flag = True
            else:
                flag = False
            fresult = DVResult(
                flag=flag,
                code=code,
                alt=alt
            )

        return sanitised_reads_by_sample, fresult
