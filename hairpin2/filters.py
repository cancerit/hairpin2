from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import ClassVar, override, Generic, TypeVar, final, cast, Any
from collections.abc import Sequence
from pysam import AlignedSegment
from hairpin2 import ref2seq as r2s
from enum import IntEnum, auto, EnumMeta
from statistics import median, stdev
# pyright: reportExplicitAny=false

# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust
# TODO: if you want to add another filter guide (commented out example minimal definition)

# TODO: method[[...], str] for returning string to write for each filter
# TODO: provide an explanation of all this strict-typing magic
# e.g. providing both static and runtime enforcment via Generic[T] and CodeEnum respectively
T = TypeVar("T", bound=IntEnum)
@dataclass
class FilterData(ABC, Generic[T]):
    CodeEnum: ClassVar[EnumMeta]
    name: str

    def __post_init__(self) -> None:
        self._flag: bool | None = None
        self._code: T | None = None

    @property
    def flag(self) -> bool | None:
        return self._flag

    @flag.setter
    def flag(self, new_state: bool):
        if not isinstance(new_state, bool):  # pyright: ignore[reportUnnecessaryIsInstance]
            raise ValueError('The flag field may only be set to a boolean value')  # pyright: ignore[reportUnreachable]
        else:
            self._flag = new_state

    @property
    def code(self) -> T | None:
        return self._code

    @code.setter
    def code(self, new_code: T) -> None:
        if not isinstance(new_code, self.CodeEnum):
            raise ValueError(f"The code field may only be filled by values given in the associated CodeEnum ClassVar")
        else:
            self._code = new_code

    @abstractmethod
    def test(
        self,
        *args,  # pyright: ignore[reportMissingParameterType, reportUnknownParameterType]
        **kwargs  # pyright: ignore[reportMissingParameterType, reportUnknownParameterType]
    ) -> None | dict[Any, list[AlignedSegment]]:
        """
        each filter must define a test method
        """


class ADFCodes(IntEnum):
    INSUFFICIENT_READS = 0
    SIXTYAI = auto()
    SIXTYBI = auto()
@final
@dataclass
class ADFilter(FilterData[ADFCodes]):
    CodeEnum: ClassVar[type[ADFCodes]] = ADFCodes
    name: str = field(default='ADF', init=False)
    variant_start: int
    edge_definition: float = 0.15  # relative proportion, by percentage, of a read to be considered 'the edge'
    edge_clustering_threshold: float = 0.9  # percentage threshold
    min_MAD_one_strand: int = 0  # exclusive (and subsequent params)
    min_sd_one_strand: float = 4
    min_MAD_both_strand_weak: int = 2
    min_sd_both_strand_weak: float = 2
    min_MAD_both_strand_strong: int = 1
    min_sd_both_strand_strong: float = 10
    min_reads: int = 1  # inclusive

    # per paper, can set hairpin for mutations distant alignment start
    # in the case where both strands have sufficient supporting reads
    @override
    def test(
        self,
        variant_start: int,
        variant_reads: Sequence[AlignedSegment],
    ) -> None:

        if len(variant_reads) < 1:
            self.CodeEnum.INSUFFICIENT_READS
            self.flag = False
            return

        # *l*engths of *a*lignment starts *to* *m*utant query positions
        la2ms_f: list[int] = []
        la2ms_r: list[int] = []
        near_start_f: list[bool] = []
        near_start_r: list[bool] = []

        for read in variant_reads:
            mut_qpos = r2s.ref2querypos(read, variant_start)
            if read.flag & 0x10:
                # +1 to include last base in length
                la2m = read.query_alignment_end - mut_qpos + 1
                near_start_r.append(((la2m / read.query_alignment_length)
                                     <= self.edge_definition))
                la2ms_r.append(la2m)
            else:
                la2m = mut_qpos - read.query_alignment_start + 1
                near_start_f.append(((la2m / read.query_alignment_length)
                                     <= self.edge_definition))
                la2ms_f.append(la2m)

        # hairpin conditions from Ellis et al. 2020, Nature Protocols
        # sometimes reported as 2021
        if len(la2ms_f) <= self.min_reads and len(la2ms_r) <= self.min_reads:
            self.code = self.CodeEnum.INSUFFICIENT_READS  # TODO: make granular
            self.flag = False
        else:
            if len(la2ms_f) > self.min_reads:  # if this, then calculate stats
                med_f = median(la2ms_f)  # range calculation replaced with true MAD calc (for r strand also)
                mad_f = median(map(lambda x: abs(x - med_f), la2ms_f))
                sd_f = stdev(la2ms_f)
                if len(la2ms_r) <= self.min_reads:  # if also this, test
                    if (((sum(near_start_f) / len(near_start_f)) < self.edge_clustering_threshold) and
                        mad_f > self.min_MAD_one_strand and
                            sd_f > self.min_sd_one_strand):
                        self.code = self.CodeEnum.SIXTYAI  # 60A(i)
                        self.flag = False
                    else:
                        self.code = self.CodeEnum.SIXTYAI
                        self.flag = True
            # the nested if statement here makes the combined condition mutually exclusive with the above
            if len(la2ms_r) > self.min_reads:
                med_r = median(la2ms_r)
                mad_r = median(map(lambda x: abs(x - med_r), la2ms_r))
                sd_r = stdev(la2ms_r)
                if len(la2ms_f) <= self.min_reads:
                    if (((sum(near_start_r) / len(near_start_r)) < self.edge_clustering_threshold) and
                        mad_r > self.min_MAD_one_strand and
                            sd_r > self.min_sd_one_strand):
                        self.code = self.CodeEnum.SIXTYAI
                        self.flag = False
                    else:
                        self.code = self.CodeEnum.SIXTYAI
                        self.flag = True
            if len(la2ms_f) > self.min_reads and len(la2ms_r) > self.min_reads:
                frac_lt_thresh = (sum(near_start_f + near_start_r)
                                  / (len(near_start_f) + len(near_start_r)))
                if (frac_lt_thresh < self.edge_clustering_threshold or
                    (mad_f > self.min_MAD_both_strand_weak and mad_r > self.min_MAD_both_strand_weak and sd_f > self.min_sd_both_strand_weak and sd_r > self.min_sd_both_strand_weak) or  # pyright: ignore[reportPossiblyUnboundVariable]
                    (mad_f > self.min_MAD_both_strand_strong and sd_f > self.min_sd_both_strand_strong) or  # pyright: ignore[reportPossiblyUnboundVariable]
                        (mad_r > self.min_MAD_both_strand_strong and sd_r > self.min_sd_both_strand_strong)):  # pyright: ignore[reportPossiblyUnboundVariable]
                    self.code = self.CodeEnum.SIXTYBI  # 60B(i)
                    self.flag = False
                else:
                    self.code = self.CodeEnum.SIXTYBI
                    self.flag = True


class ALFCodes(IntEnum):
    INSUFFICIENT_READS = 0
    ON_THRESHOLD = auto()
@final
@dataclass
class ALFilter(FilterData[ALFCodes]):
    CodeEnum: ClassVar[type[ALFCodes]] = ALFCodes
    name: str = field(default='ALF', init=False)
    _avg_as: float = field(init=False)
    al_thresh: float = 0.93

    @property
    def avg_as(self) -> float:
        return self._avg_as

    @avg_as.setter
    def avg_as(self, val: float) -> None:
        if not isinstance(val, float):
            raise ValueError(f"average AS must be represented by a float")
        else:
            self._avg_as = val

    @override
    def test(
        self,
        variant_reads: Sequence[AlignedSegment]
    ) -> None:
        if len(variant_reads) < 1:
            self.CodeEnum.INSUFFICIENT_READS
            self.flag = False
            return
        
        aln_scores: list[float] = []

        for read in variant_reads:
            try:
                aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
            except KeyError:
                pass
        if len(aln_scores) != 0:
            self.avg_as = median(aln_scores)
            self.code = self.CodeEnum.ON_THRESHOLD
            self.flag = False
            if self.avg_as <= self.al_thresh:
                self.flag = True
        else:
            self.code = self.CodeEnum.INSUFFICIENT_READS
            self.flag = False


class DVFCodes(IntEnum):
    INSUFFICIENT_READS = 0
    DUPLICATION = auto()
@final
@dataclass
class DVFilter(FilterData[DVFCodes]):
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
    CodeEnum: ClassVar[type[DVFCodes]] = DVFCodes
    name: str = field(default='DVF', init=False)
    min_boundary_deviation: int = 6  # TODO: can be used to turn off
    # TODO: n.b. neither of these options prevent read removal due to duplication, so test still functions as QC (and I think this is fine, just document more)
    read_number_difference_threshhold: int = 0 # change in reads! TODO: express as fraction? discuss with Peter/Phuong
    nsamples_threshold: int = 1

    # detect PCR duplicates previously missed due to slippage
    # this implementation assumes that sorting on first element of each sublist
    # is appropriate, per Peter's initial implementation.
    # is an all against all comparison between all read lists more appropriate?
    # and between pairs of readlists, why is comparing sorted pairs most appropriate?
    @override
    def test(
        self,
        variant_reads_by_sample: dict[str, list[AlignedSegment]]
    ) -> dict[str, list[AlignedSegment]]:
        """
        A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.

        """
        nsamples_with_duplication = 0
        if not any([len(reads) > 1 for reads in variant_reads_by_sample.values()]):
            self.code = self.CodeEnum.INSUFFICIENT_READS
            self.flag = False
            return variant_reads_by_sample
        else:
            self.code = self.CodeEnum.DUPLICATION  # testing possible, and this is the only code
            sanitised_reads_by_sample: dict[str, list[AlignedSegment]] = {}
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
                        if all([x <= self.min_boundary_deviation for x in max_diff_per_comparison]):
                            # then the read pair being examined is a duplicate of the others in the pool
                            dup_endpoint_test_pool.append(testing_endpoints)
                            dup_idcs.append(original_index)  # store original index of this duplicate
                        else:
                            # read at i is not dup of reads in dup_endpoint_test_pool
                            # start again, test read at i against reads subsequent from i in ends_sorted
                            dup_endpoint_test_pool = [testing_endpoints]
                    sanitised_reads = [read for i, read in enumerate(reads) if i not in dup_idcs]
                    if (len(sanitised_reads) - len(reads)) > self.read_number_difference_threshhold:
                        nsamples_with_duplication += 1
                    sanitised_reads_by_sample[sample_key] = sanitised_reads

            if nsamples_with_duplication > self.nsamples_threshold:
                self.flag = True

            return sanitised_reads_by_sample



# TODO: QCfilter after discussion with Peter/Phuong
# @dataclass
# class QCFilter(FilterData):
#     """
#     All reads supporting the variant have failed hairpin2 QC
#     """
#     name: str = field(default='QCF')


# useful aliases
type AnyFilterSequence = Sequence[FilterData[Any]]
