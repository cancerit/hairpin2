from array import array
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Flag, auto
from statistics import mean, median, stdev
from typing import Any, cast, final

from pysam import AlignedSegment

from hairpin2.const import MutTypes, Strand, TaggerNamespaces, Tags
from hairpin2.infrastructure.structures import (
    ExtendedRead,
    TestOutcomes,
    mark_read,
    read_has_mark,
    record_operation,
)
from hairpin2.utils.ref2seq import CigarError, ref2querypos, ref_end_via_cigar

# DEV: pysam has issues with its type hints, so:
# pyright: reportUnknownVariableType = false
# pyright: reportUnknownMemberType = false
# pyright: reportUnknownArgumentType = false
# TODO (DEV):
# define a protocol for taggers


# TODO:
# Make a global filter for those properties which all other processes need validated to function, excluding reads without those properties

"""
Classes and functions for additive read tagging processes that apply various tags to reads for subsequent filtering and analysis.
"""


def _sort_reads_by_qual[T: ExtendedRead | AlignedSegment](reads: Iterable[T]) -> list[T]:
    """
    Internal function for use when needing to select the best read of a pool (2 or more reads), to ensure this is always done in the same way.
    """
    return sorted(reads, key=lambda x: mean(x.query_qualities))  # pyright: ignore[reportArgumentType]


class TagSupportingReads:
    """
    pc8's Support Assessor

    Confirm whether a read supports a variant and optionally mark supporting reads.

    Used to tag reads as supporting a variant. This tag can then be used by other processes
    to include/exclude reads with that tag from their analysis.
    """

    class SupportFlags(Flag):
        CLEAR = 0
        FIELDS_MISSING = auto()
        NOT_ALIGNED = auto()
        NOT_ALT = auto()
        SHORT = auto()
        BAD_OP = auto()

    @classmethod
    def check_read_supporting(
        cls,
        read: ExtendedRead,
        mut_type: MutTypes,
        alt: str,
        vcf_start: int,
        vcf_stop: int,
        mark: bool = True,
    ) -> bool:
        if mut_type not in MutTypes:
            raise ValueError(
                f"Unsupported mutation type {mut_type} provided to mut_type argument - supports {MutTypes} only"
            )

        # if this check is passed, the subsequent
        # two conditions can only be error states
        # hence raise rather than return false
        if not ((read.flag & 3) == 3) or read.flag & 3852:
            return False

        if any(
            flg is None
            for flg in [
                read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
            ]
        ) or not read.has_tag("MC"):
            raise ValueError(f"Essential fields missing from alignment record {read.query_name}, data in alignment file corrupted")

        tag_dat = read.get_tag("MC", with_value_type=True)
        if (tag_dat[1] != 'Z'):
            raise RuntimeError(f"{read.query_name} MC tag is not of type 'Z' as mandated by SAM format spec. Data in alignment file corrupted.")

        support_flag = cls.SupportFlags.CLEAR
        if mut_type in [MutTypes.SUB, MutTypes.INS]:
            try:
                mut_pos = ref2querypos(read, vcf_start)
            except ValueError:
                support_flag |= cls.SupportFlags.NOT_ALIGNED
            else:
                if cast(str, read.query_sequence)[mut_pos : mut_pos + len(alt)] != alt:
                    support_flag |= cls.SupportFlags.NOT_ALT
                if (
                    mut_type == MutTypes.INS
                ):  # INS - mut_pos is position immediately before insertion
                    if mut_pos + len(alt) > read.query_length:
                        support_flag |= cls.SupportFlags.SHORT
                    else:
                        mut_alns = [
                            (q, r)
                            for q, r in read.get_aligned_pairs()
                            if q in range(mut_pos + 1, mut_pos + len(alt))
                        ]
                        if any(r is not None for _, r in mut_alns):
                            support_flag |= cls.SupportFlags.BAD_OP
        elif mut_type == MutTypes.DEL:
            rng = list(range(vcf_start, vcf_stop + 1))
            mut_alns = [q for q, r in read.get_aligned_pairs() if r in rng]
            if len(mut_alns) != len(rng):
                support_flag |= cls.SupportFlags.SHORT
            if any(x is not None for x in mut_alns[1:-1]) or any(
                x is None for x in [mut_alns[0], mut_alns[-1]]
            ):
                support_flag |= cls.SupportFlags.BAD_OP

        # ValidatorFlag not returned but kept in case useful in the future
        if support_flag == cls.SupportFlags.CLEAR:
            support = True
            if mark:
                mark_read(read, Tags.SUPPORT_TAG)
        else:
            support = False

        record_operation(read, TaggerNamespaces.MARK_SUPPORT)
        return support


class TagLowQualReads:
    """
    pc8's Quality Assessor

    Determine if a read is of low quality and optionally mark low-quality reads.

    The result of this process is used to mark reads as low quality. This tag can then be used by other processes
    to include/exclude reads with that tag from their analysis. The fraction of reads tagged as low quality,
    or as duplicates by the duplicate marking function below, is used to call the LQF flag on variants
    """

    class QualFlags(Flag):
        CLEAR = 0
        MAPQUAL = auto()
        CLIPQUAL = auto()
        BASEQUAL = auto()

    @classmethod
    def check_low_qual_read(
        cls,
        read: ExtendedRead,
        vcf_start: int,
        alt: str,
        mut_type: MutTypes,
        min_basequal: int,
        min_mapqual: int,
        min_clipqual: int,
        mark: bool = True,
    ) -> bool:
        qual_flag = cls.QualFlags.CLEAR  # 0 - evaluates false

        if read.mapping_quality < min_mapqual:
            qual_flag |= cls.QualFlags.MAPQUAL

        if (
            "S" in read.cigarstring  # pyright: ignore[reportOperatorIssue]
            and mean(read.query_alignment_qualities) < min_clipqual  # pyright: ignore[reportArgumentType]
        ):
            qual_flag |= cls.QualFlags.CLIPQUAL

        if mut_type == MutTypes.SUB:
            try:
                mut_pos = ref2querypos(read, vcf_start)
            except ValueError:
                raise ValueError(
                    f"read: {read.query_name} appears not to cover variant at {vcf_start} - this function should only be run on reads confirmed to be supporting the variant"
                )
            else:
                if any(
                    bq < min_basequal
                    for bq in cast(array[Any], read.query_qualities)[mut_pos : mut_pos + len(alt)]
                ):
                    qual_flag |= cls.QualFlags.BASEQUAL

        if qual_flag == cls.QualFlags.CLEAR:
            low_qual = False
        else:
            low_qual = True
            if mark:
                mark_read(read, Tags.LOW_QUAL_TAG)

        record_operation(read, TaggerNamespaces.MARK_LOW_QUAL)
        return low_qual


class TagFragmentReads:
    """
    ab63 Overlap Assessor

    Mark one read from a fragment pair if both members of the fragment pair are present in the input reads.

    Used to tag reads as an overlapping mate. This tag can then be used by other processes to
    include/exclude reads with that tag from their analysis
    """

    @staticmethod
    def check_for_mates(reads: Iterable[ExtendedRead]):
        """
        Mark lower mean qual member of a read pair as overlapping
        """
        seen_qnames: dict[str, ExtendedRead] = {}
        for read in reads:
            qn = read.query_name
            rq = read.query_qualities
            if qn is None or rq is None:
                continue
            if qn in seen_qnames:
                sorted_reads = _sort_reads_by_qual([read, seen_qnames[qn]])
                mark_read(sorted_reads[0], Tags.OVERLAP_TAG)  # mark lower qual
            else:
                seen_qnames[qn] = read

            record_operation(read, TaggerNamespaces.MARK_OVERLAP)


class TagStutterDuplicateReads:
    """
    pc8's Stutter Duplicates Assessor

    Tag hidden duplicate reads which have shifted endpoints due to PCR stutter (hence evading normal dupmarking).

    The result of this process is used to mark reads as duplicates. This tag can then be used by other processes
    to include/exclude reads with that tag from their analysis. The fraction of reads tagged as duplicates
    by this process is used to assess a variant for the DVF flag and, partly (with the low-quality tag), the LQF flag
    """

    @dataclass
    class FragmentEndpoints:
        read: ExtendedRead
        endpoints: tuple[int, int, int, int]

    @classmethod
    def check_stutter_duplicates(
        cls, reads: Iterable[ExtendedRead], duplication_window_size: int = 6
    ):
        if not reads:
            return

        # prep data
        endpointsl = []
        for read in reads:
            try:
                tag_dat = read.get_tag("MC", with_value_type=True)
                if (tag_dat[1] != 'Z'):
                    raise RuntimeError("MC tag is not of type 'Z' as mandated by SAM format spec. Data in alignment file corrupted.")
                mate_cig = str(tag_dat[0])
                next_ref_end = ref_end_via_cigar(mate_cig, read.next_reference_start)
            except (KeyError, CigarError):
                continue
            if read.reference_end is None:
                continue
            if read.query_qualities is None:
                continue

            pair_endpoints: tuple[int, int, int, int] = cast(
                tuple[int, int, int, int],
                tuple(
                    sorted(
                        [
                            read.reference_start,
                            read.reference_end,
                            read.next_reference_start,
                            next_ref_end,
                        ]
                    )
                ),
            )
            endpointsl.append(cls.FragmentEndpoints(read, pair_endpoints))
        # sort by starting endpoint
        endpointsl = sorted(endpointsl, key=lambda x: x.endpoints[0])

        # test data
        # NOTE: do you want to only look at supporting reads for this dupmarking? Discuss
        # setup test
        init_case = endpointsl.pop(0)
        dup_endpoint_comparison_pool: list[tuple[int, int, int, int]] = []
        dup_endpoint_comparison_pool.append(init_case.endpoints)
        dup_read_pool: list[ExtendedRead] = [init_case.read]
        record_operation(init_case.read, TaggerNamespaces.MARK_STUTTER_DUP)

        for i in range(len(endpointsl)):
            read = endpointsl[i].read
            testing_endpoints = endpointsl[i].endpoints

            max_diff_per_comparison: list[int] = []
            for comparison_endpoints in dup_endpoint_comparison_pool:
                endpoint_diffs: list[int] = [
                    abs(x - y) for x, y in zip(comparison_endpoints, testing_endpoints)
                ]
                # assert len(endpoint_diffs) == 4  # sanity check - commented out but left in as a hint
                max_diff_per_comparison.append(max(endpoint_diffs))

            if all(x <= duplication_window_size for x in max_diff_per_comparison):
                # then the read pair being examined is a duplicate of the others in the pool
                dup_endpoint_comparison_pool.append(
                    testing_endpoints
                )  # add the endpoints to the comparison pool
                dup_read_pool.append(read)  # add the read to the dup list
            else:
                assert len(dup_endpoint_comparison_pool) == len(dup_read_pool)

                # read at i is not dup of reads in dup_endpoint_test_pool
                # since reads are sorted by endpoints, no more dups incoming
                # so finalise the current dup pool
                if len(dup_read_pool) > 1:
                    sorted_dups = _sort_reads_by_qual(dup_read_pool)
                    _ = sorted_dups.pop()  # pop the highest quality
                    for read in sorted_dups:
                        mark_read(read, Tags.STUTTER_DUP_TAG)  # mark the others

                # start again, test the read at i against reads subsequent from i in ends_sorted
                # reset the pools
                dup_read_pool = [read]
                dup_endpoint_comparison_pool = [testing_endpoints]
            record_operation(read, TaggerNamespaces.MARK_STUTTER_DUP)


"""
Classes and functions for read-aware variant flagging tests.
"""


@final
class AlignmentScoreTest:
    """
    Pass/Fail a variant on the average alignment score of supporting reads.
    """

    @dataclass
    class ResultPack:
        class Info(Flag):
            NODATA = 0
            NO_READS = 1
            INSUFFICIENT_AS_TAGS = 2
            ON_THRESHOLD = 4

        outcome: TestOutcomes
        avg_as: float | None
        reason: Info

    # Alignment Score Test -- ALF flag
    # flags must be powers of 2

    @classmethod
    def test_variant_reads(
        cls, reads: Iterable[ExtendedRead], avg_AS_threshold: float
    ) -> ResultPack:
        """
        pc8's Alignment Score Assessor
        """
        if not reads:
            result = cls.ResultPack(TestOutcomes.NA, None, cls.ResultPack.Info.NO_READS)
        else:
            aln_scores: list[float] = []

            for read in reads:
                try:
                    tag_dat = read.get_tag("AS", with_value_type=True)
                    if (tag_dat[1] not in ['i', 'C']): # BUG/TODO revisit mysterious undescribed "C" type, appears to be an idiosyncracy of htslib at the leaset
                        breakpoint()
                        raise RuntimeError(f"{read.query_name} AS tag is not of type 'i' as mandated by SAM format spec. Data in alignment file corrupted.")
                    ascore = int(tag_dat[0])
                    aln_scores.append(ascore / read.query_length)
                except KeyError:
                    pass
            if len(aln_scores) != 0:
                avg_as = median(aln_scores)
                if avg_as <= avg_AS_threshold:
                    outcome = TestOutcomes.VARIANT_FAIL
                else:
                    outcome = TestOutcomes.VARIANT_PASS

                result = cls.ResultPack(outcome, avg_as, cls.ResultPack.Info.ON_THRESHOLD)
            else:
                result = cls.ResultPack(
                    TestOutcomes.NA, None, cls.ResultPack.Info.INSUFFICIENT_AS_TAGS
                )

        return result


# Anomalous Distribution Test -- ADF flag
@final
class AnomalousDistributionTest:
    """
    Pass/Fail a variant on the distribution of the mutant base position on the supporting reads

    Uses the conditions given by Ellis et al. in `doi.org/10.1038/s41596-020-00437-6 <https://doi.org/10.1038/s41596-020-00437-6>`_,
    and one additional condition devised by al35

    The full text of the conditional is reproduced below, with editorials in []. There is a point of ambiguity in the original conditional.
    The interpretation that this tool has opted for is indicated by [] and is expanded upon subsequently.

    ---------------------------------------------------------

    For each variant, if the number of variant-supporting reads determined [IN PRIOR STEPS] is low (i.e. 0–1
    reads) for one strand, follow Option A. For each variant, if both strands have sufficient variant-supported reads
    (i.e. ≥2 reads), follow Option B.

    A) Low number of variant-supporting reads on one strand
        (i) For each variant, if one strand had too few variant-supporting reads, the other strand must
        conform to:
            - Fewer than 90% of variant-supporting reads [ON THE STRAND] have the variant located within the first 15% of the read measured from the alignment start position.
            - MAD >0 and s.d. >4 for that strand.

    B) Sufficient variant-supporting reads on both strands
        (i) For each variant, if both strands have sufficient variant-supporting reads (i.e., ≥2 reads),
        then one of the following must be true:
            - Fewer than 90% of variant-supporting reads should have the variant located within the first 15% of the read measured from the alignment start position.
            - MAD >2 and s.d. >2 for both strands.
            - MAD >1 and s.d. >10 for one strand (i.e., strong evidence of high variability in variant position in variant-supporting reads).
    ---------------------------------------------------------

    The point of ambiguity is whether, on path A, to include the single read from the strand which does not sufficiently
    support the variant in the test of positional distribution across the supporting reads.
    The present interpretation is that the test should be applied only to the reads on the "other strand",
    since the phrasing "the other strand must conform to" implies the exclusion of the single read on the low-support
    strand.

    The additional condition to be checked is simply whether at least N supporting reads express away from the read
    edge. This condition is experimental at this time and can be disabled by setting min_non_edge_reads to 0
    """

    @dataclass
    class ResultPack:
        class Info(Flag):
            NODATA = 0
            NO_READS = 1
            INSUFFICIENT_READS = 2
            EDGE_CLUSTERING = 4
            ONE_STRAND_DISTRIB = 8
            BOTH_STRAND_DISTRIB_BOTH = 16
            BOTH_STRAND_DISTRIB_ONE = 32
            MIN_NON_EDGE = 64

        outcome: TestOutcomes
        strand: Strand
        reason: Info

    # NOTE:
    # per paper, can set hairpin for mutations distant from alignment start
    # in the case where both strands have sufficient supporting reads
    @classmethod
    def test_variant_reads(
        cls,
        reads: Iterable[ExtendedRead],
        record_start: int,  # VCF record start
        edge_definition: float,  # relative proportion, by percentage, of a read to be considered 'the edge'
        edge_clustering_threshold: float,  # percentage threshold
        min_MAD_one_strand: float,  # exclusive (and subsequent params)
        min_sd_one_strand: float,
        min_MAD_both_strand_weak: float,
        min_sd_both_strand_weak: float,
        min_MAD_both_strand_strong: float,
        min_sd_both_strand_strong: float,
        low_n_supporting_reads_boundary: int,  # inclusive
        min_non_edge_reads: int,
    ) -> ResultPack:
        if not reads:
            result = cls.ResultPack(
                TestOutcomes.NA,
                Strand.BOTH,
                cls.ResultPack.Info.NO_READS,
            )
        else:
            # *l*engths of *a*lignment starts *to* *m*utant query positions
            la2ms_f: list[int] = []
            la2ms_r: list[int] = []
            near_start_f: list[bool] = []
            near_start_r: list[bool] = []

            for read in reads:
                try:
                    mut_qpos = ref2querypos(
                        read, record_start
                    )  # TODO: this is additive processing and should be eventually separated as such
                except ValueError:
                    raise ValueError(f"read {read.query_name} does not cover variant")

                if read.flag & 0x10:
                    # +1 to include the last base in length
                    la2m = read.query_alignment_end - mut_qpos + 1
                    near_start_r.append((la2m / read.query_alignment_length) <= edge_definition)
                    la2ms_r.append(la2m)
                else:
                    la2m = mut_qpos - read.query_alignment_start + 1
                    near_start_f.append((la2m / read.query_alignment_length) <= edge_definition)
                    la2ms_f.append(la2m)

            # hairpin conditions from Ellis et al. 2020, Nature Protocols
            # sometimes reported as 2021
            if (
                len(la2ms_f) <= low_n_supporting_reads_boundary
                and len(la2ms_r) <= low_n_supporting_reads_boundary
            ):
                result = cls.ResultPack(
                    TestOutcomes.NA, Strand.BOTH, cls.ResultPack.Info.INSUFFICIENT_READS
                )
            else:
                strand = None
                outcome = TestOutcomes.VARIANT_PASS
                info_bits = cls.ResultPack.Info.NODATA
                if len(la2ms_f) > low_n_supporting_reads_boundary:  # if true, calculate stats
                    med_f = median(
                        la2ms_f
                    )  # range calculation replaced with true MAD calc (for r strand also)
                    mad_f = median(abs(val - med_f) for val in la2ms_f)
                    sd_f = stdev(la2ms_f)
                    if (
                        len(la2ms_r) <= low_n_supporting_reads_boundary
                    ):  # if also this, test on this path
                        strand = Strand.F
                        cond_edge_clustering = (
                            sum(near_start_f) / len(near_start_f)
                        ) < edge_clustering_threshold
                        cond_one_strand_mad = mad_f > min_MAD_one_strand
                        cond_one_strand_sd = sd_f > min_sd_one_strand

                        # this could be more concise, but I think this is the most readable approach I've found
                        if not cond_edge_clustering:
                            info_bits |= cls.ResultPack.Info.EDGE_CLUSTERING
                            outcome = TestOutcomes.VARIANT_FAIL
                        if not (cond_one_strand_mad and cond_one_strand_sd):
                            info_bits |= cls.ResultPack.Info.ONE_STRAND_DISTRIB
                            outcome = TestOutcomes.VARIANT_FAIL

                        if outcome == TestOutcomes.VARIANT_PASS:
                            info_bits = (
                                cls.ResultPack.Info.ONE_STRAND_DISTRIB
                                | cls.ResultPack.Info.EDGE_CLUSTERING
                            )
                # the nested if statement here makes the combined condition mutually exclusive with the above
                if len(la2ms_r) > low_n_supporting_reads_boundary:
                    med_r = median(la2ms_r)
                    mad_r = median(abs(val - med_r) for val in la2ms_r)
                    sd_r = stdev(la2ms_r)
                    if len(la2ms_f) <= low_n_supporting_reads_boundary:
                        strand = Strand.R
                        outcome = TestOutcomes.VARIANT_PASS
                        info_bits = cls.ResultPack.Info.NODATA
                        cond_edge_clustering = (
                            sum(near_start_r) / len(near_start_r)
                        ) < edge_clustering_threshold
                        cond_one_strand_mad = mad_r > min_MAD_one_strand
                        cond_one_strand_sd = sd_r > min_sd_one_strand

                        if not cond_edge_clustering:
                            info_bits |= cls.ResultPack.Info.EDGE_CLUSTERING
                            outcome = TestOutcomes.VARIANT_FAIL
                        if not (cond_one_strand_mad and cond_one_strand_sd):
                            info_bits |= cls.ResultPack.Info.ONE_STRAND_DISTRIB
                            outcome = TestOutcomes.VARIANT_FAIL

                        if outcome == TestOutcomes.VARIANT_PASS:
                            info_bits = (
                                cls.ResultPack.Info.ONE_STRAND_DISTRIB
                                | cls.ResultPack.Info.EDGE_CLUSTERING
                            )
                if (
                    len(la2ms_f) > low_n_supporting_reads_boundary
                    and len(la2ms_r) > low_n_supporting_reads_boundary
                ):
                    strand = Strand.BOTH
                    outcome = TestOutcomes.VARIANT_PASS
                    info_bits = cls.ResultPack.Info.NODATA
                    frac_lt_thresh = sum(near_start_f + near_start_r) / (
                        len(near_start_f) + len(near_start_r)
                    )

                    cond_edge_clustering = frac_lt_thresh < edge_clustering_threshold
                    cond_both_strand_distrib_both = (
                        mad_f > min_MAD_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                        and mad_r > min_MAD_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                        and sd_f > min_sd_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                        and sd_r > min_sd_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                    )
                    cond_both_strand_distrib_one = (
                        (mad_f > min_MAD_both_strand_strong and sd_f > min_sd_both_strand_strong)  # pyright: ignore[reportPossiblyUnboundVariable]
                        or (mad_r > min_MAD_both_strand_strong and sd_r > min_sd_both_strand_strong)  # pyright: ignore[reportPossiblyUnboundVariable]
                    )

                    if cond_edge_clustering:  # if satisfied, add cond
                        info_bits |= cls.ResultPack.Info.EDGE_CLUSTERING
                    if cond_both_strand_distrib_both:
                        info_bits |= cls.ResultPack.Info.BOTH_STRAND_DISTRIB_BOTH
                    if cond_both_strand_distrib_one:
                        info_bits |= cls.ResultPack.Info.BOTH_STRAND_DISTRIB_ONE

                    # if no pass conditions are satisfied
                    if info_bits == cls.ResultPack.Info.NODATA:
                        outcome = TestOutcomes.VARIANT_FAIL
                        info_bits |= (
                            cls.ResultPack.Info.EDGE_CLUSTERING
                            | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_BOTH
                            | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_ONE
                        )

                reads_not_near_edge = (len(near_start_f) - sum(near_start_f)) + (
                    len(near_start_r) - sum(near_start_r)
                )
                if not reads_not_near_edge > min_non_edge_reads:
                    if outcome == TestOutcomes.VARIANT_PASS:
                        outcome = TestOutcomes.VARIANT_FAIL
                        info_bits = cls.ResultPack.Info.MIN_NON_EDGE  # we failed here
                    else:
                        info_bits |= cls.ResultPack.Info.MIN_NON_EDGE  # add to other failures

                # add implied passes if total success
                if outcome == TestOutcomes.VARIANT_PASS:
                    info_bits |= (
                        cls.ResultPack.Info.MIN_NON_EDGE
                        | cls.ResultPack.Info.INSUFFICIENT_READS
                        | cls.ResultPack.Info.NO_READS
                    )

                assert strand is not None  # appease a type checker

                result = cls.ResultPack(outcome, strand, info_bits)

        return result


class ProportionBasedTest:
    """
    Pass/Fail a variant on a proportion of supporting reads with/without a given property.

    Used by the LQF flag test to test a proportion of stutter duplicate and low-quality supporting reads.
    Used by the DVF flag test to test a proportion of stutter duplicate reads.
    """

    @dataclass
    class ResultPack:
        class Info(Flag):
            NODATA = 0
            NO_READS = 1
            THRESHOLD = 2
            MIN_PASS = 4

        outcome: TestOutcomes
        reason: Info
        prop_loss: float

    @classmethod
    def test_variant_reads(
        cls,
        reads: Sequence[ExtendedRead],
        tags_to_check: Iterable[Tags],
        read_loss_threshold: float,
        min_without: int = 0,
    ) -> ResultPack:
        if not reads:
            result = cls.ResultPack(TestOutcomes.NA, cls.ResultPack.Info.NO_READS, 0)
        else:
            info_bits = cls.ResultPack.Info.NODATA
            outcome = TestOutcomes.VARIANT_PASS
            ntotal = len(reads)

            nwith = 0
            for read in reads:
                nwith += int(bool(sum(read_has_mark(read, tag) for tag in tags_to_check)))
            nwithout = ntotal - nwith
            loss_ratio = (
                nwith / ntotal
            )  # you could alternately formulate this as nwith / nwithout...
            if loss_ratio > read_loss_threshold:
                info_bits |= cls.ResultPack.Info.THRESHOLD
                outcome = TestOutcomes.VARIANT_FAIL
            if nwithout < min_without:
                info_bits |= cls.ResultPack.Info.MIN_PASS
                outcome = TestOutcomes.VARIANT_FAIL

            if outcome == TestOutcomes.VARIANT_PASS:
                info_bits = ~cls.ResultPack.Info(0)  # satisfied all conditions

            result = cls.ResultPack(outcome, info_bits, loss_ratio)

        return result
