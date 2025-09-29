from array import array
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from enum import Flag
from statistics import mean, median, stdev
from typing import Any, cast, final

from hairpin2.const import MutTypes, Strand, Tags, ValidatorFlags
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


# SECTION --- functions/classes used by additive read tagging processes to add various tags to reads for subsequent filtering and analysis.


# Used to tag reads as supporting a variant
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis
class TagSupportingReads:
    """
    pc8's Support Assessor

    Confirm whether a read supports a variant, and optionally mark supporting reads.
    """

    @staticmethod
    def check_read_supporting(
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

        support_flag = ValidatorFlags.CLEAR
        try:
            mate_cig = str(read.get_tag("MC"))
        except KeyError:
            mate_cig = None
        else:
            if (
                not mate_cig[0].isdigit()
                or not all([(c.isalnum() or c == "=") for c in mate_cig])
                or len(mate_cig) < 2
            ):
                mate_cig = None
        if any(
            flg is None
            for flg in [
                read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig,
            ]
        ):
            support_flag |= ValidatorFlags.READ_FIELDS_MISSING
        elif mut_type in [MutTypes.SUB, MutTypes.INS]:
            try:
                mut_pos = ref2querypos(read, vcf_start)
            except ValueError:
                support_flag |= ValidatorFlags.NOT_ALIGNED
            else:
                if mut_type == MutTypes.SUB:
                    if cast(str, read.query_sequence)[mut_pos : mut_pos + len(alt)] != alt:
                        support_flag |= ValidatorFlags.NOT_ALT
                if (
                    mut_type == MutTypes.INS
                ):  # INS - mut_pos is position immediately before insertion
                    if mut_pos + len(alt) > read.query_length:
                        support_flag |= ValidatorFlags.SHORT
                    else:
                        mut_alns = [
                            (q, r)
                            for q, r in read.get_aligned_pairs()
                            if q in range(mut_pos + 1, mut_pos + len(alt) + 1)
                        ]
                        if any([r is not None for _, r in mut_alns]):
                            support_flag |= ValidatorFlags.BAD_OP
                        if (
                            cast(str, read.query_sequence)[mut_pos + 1 : mut_pos + len(alt) + 1]
                            != alt
                        ):
                            support_flag |= ValidatorFlags.NOT_ALT
        elif mut_type == MutTypes.DEL:
            rng = list(range(vcf_start, vcf_stop + 1))
            mut_alns = [q for q, r in read.get_aligned_pairs() if r in rng]
            if len(mut_alns) != len(rng):
                support_flag |= ValidatorFlags.SHORT
            if any([x is not None for x in mut_alns[1:-1]]) or any(
                [x is None for x in [mut_alns[0], mut_alns[-1]]]
            ):
                support_flag |= ValidatorFlags.BAD_OP

        # ValidatorFlag not returned, but kept in case useful in future
        if support_flag == ValidatorFlags.CLEAR:
            support = True
            if mark:
                mark_read(read, Tags.SUPPORT_TAG)
        else:
            support = False

        record_operation(read, "mark-support")
        return support


# The result of this process is used to mark reads as low quality.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis.
# The fraction of reads tagged as low quality, or as duplicates
# by the duplicate marking function below, is used to call the
# LQF flag on variants
class TagLowQualReads:
    """
    pc8's Quality Assessor

    Determine if a read is of low quality, and optionally mark low quality reads.
    """

    @staticmethod
    def check_low_qual_read(
        read: ExtendedRead,
        vcf_start: int,
        alt: str,
        mut_type: MutTypes,
        min_basequal: int,
        min_mapqual: int,
        min_clipqual: int,
        mark: bool = True,
    ) -> bool:
        qual_flag = ValidatorFlags.CLEAR  # 0 - evaluates false

        if not (read.flag & 0x2) or read.flag & 0xE00:
            qual_flag |= ValidatorFlags.FLAG

        if read.mapping_quality < min_mapqual:
            qual_flag |= ValidatorFlags.MAPQUAL

        if (
            "S" in read.cigarstring  # pyright: ignore[reportOperatorIssue]
            and mean(read.query_alignment_qualities) < min_clipqual  # pyright: ignore[reportArgumentType]
        ):
            qual_flag |= ValidatorFlags.CLIPQUAL

        if mut_type == MutTypes.SUB:
            try:
                mut_pos = ref2querypos(read, vcf_start)
            except ValueError:
                raise ValueError(
                    f"read: {read.query_name} appears not to cover variant at {vcf_start} - this function should only be run on reads confirmed to be supporting the variant"
                )
            else:
                if any(
                    [
                        bq < min_basequal
                        for bq in cast(array[Any], read.query_qualities)[
                            mut_pos : mut_pos + len(alt)
                        ]
                    ]
                ):
                    qual_flag |= ValidatorFlags.BASEQUAL

        if qual_flag == ValidatorFlags.CLEAR:
            low_qual = False
        else:
            low_qual = True
            if mark:
                mark_read(read, Tags.LOW_QUAL_TAG)

        record_operation(read, "mark-low-qual")
        return low_qual


# Used to tag reads as an overlapping mate.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis
class TagFragmentOverlapReads:
    """
    pc8's Overlap Assessor

    Determine if a read is an overlapping mate for this position, and optionally mark overlapping mates.
    """

    @staticmethod
    def check_fragment_overlap(read: ExtendedRead, vcf_start: int, mark: bool = True):
        overlap = False

        mate_cig = str(read.get_tag("MC"))  # will error if no tag

        # NOTE: introduces strand bias!!
        # NOTE: does not check if overlapping member supports variant!
        if read.flag & 0x80:  # if second in pair
            read_range = range(
                read.reference_start,
                read.reference_end,  # pyright: ignore[reportArgumentType]
            )
            mate_range = range(
                read.next_reference_start, ref_end_via_cigar(mate_cig, read.next_reference_start)
            )
            overlapping_positions = set(read_range).intersection(mate_range)
            if vcf_start in overlapping_positions:
                overlap = True

        if overlap and mark:
            mark_read(read, Tags.OVERLAP_TAG)

        record_operation(read, "mark-overlap")

        return overlap


# The result of this process is used to mark reads as duplicates.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis.
# The fraction of reads tagged as duplicates by this process is
# used to assess a variant for the DVF flag
# and, partly (with the low quality tag), the LQF flag
class TagStutterDuplicateReads:
    """
    pc8's Stutter Duplicates Assessor

    Tag hidden duplicate reads which have shifted endpoints due to PCR stutter (hence evading normal dupmarking).
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
                mate_cig = read.get_tag("MC")
                next_ref_end = ref_end_via_cigar(str(mate_cig), read.next_reference_start)
            except (KeyError, CigarError):
                continue
            if read.reference_end is None:
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
        endpointsl = sorted(endpointsl, key=lambda x: x.endpoints[0])

        # test data
        # NOTE: do you want to only look at supporting reads for this dupmarking? Discuss
        # setup test
        init_case = endpointsl.pop(0)
        dup_endpoint_comparison_pool: list[tuple[int, int, int, int]] = []
        dup_endpoint_comparison_pool.append(init_case.endpoints)
        dup_read_pool: list[ExtendedRead] = [init_case.read]
        record_operation(init_case.read, "mark-duplicates")

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

            if all([x <= duplication_window_size for x in max_diff_per_comparison]):
                # then the read pair being examined is a duplicate of the others in the pool
                dup_endpoint_comparison_pool.append(
                    testing_endpoints
                )  # add the endpoints to the comparison pool
                dup_read_pool.append(read)  # add the read to the dup list
            else:
                assert len(dup_endpoint_comparison_pool) == len(dup_read_pool)

                # read at i is not dup of reads in dup_endpoint_test_pool
                # start again, test read at i against reads subsequent from i in ends_sorted
                # NOTE: this means reads found after are not tested against previous duplicate groups...
                if len(dup_read_pool) > 1:
                    dup_read_pool.sort(
                        key=lambda x: mean(x.query_qualities)
                    )  # type checker believes there is an error here. I have not silenced it because in principle it is correct, however it is pysam's responsibility as it is a failing in their library typing.
                    _ = dup_read_pool.pop()  # pop the highest quality
                    for read in dup_read_pool:
                        mark_read(read, Tags.STUTTER_DUP_TAG)  # mark the others

                # reset the pools
                dup_read_pool = [read]
                dup_endpoint_comparison_pool = [testing_endpoints]
            record_operation(read, "mark-duplicates")


# END SECTION --- read tagging functions/classes


# SECTION --- read-aware variant flagging tests


@final
class AlignmentScoreTest:
    """
    Pass/Fail a variant on the average alignment score of supporting reads.
    """

    @dataclass
    class ResultPack:
        class Info(Flag):
            NODATA = 0
            INSUFFICIENT_READS = 1
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
            result = cls.ResultPack(TestOutcomes.NA, None, cls.ResultPack.Info.INSUFFICIENT_READS)
        else:
            aln_scores: list[float] = []

            for read in reads:
                try:
                    aln_scores.append(int(read.get_tag("AS")) / read.query_length)
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

    Uses the conditions given by Ellis et al. in doi.org/10.1038/s41596-020-00437-6, and one additional condition devised by al35

    The full text of the conditional is reproduced below, with editorials in []. There is a point of ambiguity in the original conditional.
    The interpretation that this tool has opted for is indicated by [], and is expanded upon subsequently.

    ---------------------------------------------------------

    For each variant, if the number of variant-supporting reads determined [IN PRIOR STEPS] is low (i.e., 0–1
    reads) for one strand, follow Option A. For each variant, if both strands have sufficient variant-
    supported reads (i.e., ≥2 reads), follow Option B.

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
            NO_TESTABLE_READS = 1
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
        min_MAD_one_strand: int,  # exclusive (and subsequent params)
        min_sd_one_strand: float,
        min_MAD_both_strand_weak: int,
        min_sd_both_strand_weak: float,
        min_MAD_both_strand_strong: int,
        min_sd_both_strand_strong: float,
        low_n_supporting_reads_boundary: int,  # inclusive
        min_non_edge_reads: int,
    ) -> ResultPack:
        if not reads:
            result = cls.ResultPack(
                TestOutcomes.NA,
                Strand.BOTH,
                cls.ResultPack.Info.NO_TESTABLE_READS,
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
                    # +1 to include last base in length
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
                    mad_f = median(map(lambda x: abs(x - med_f), la2ms_f))
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
                                | cls.ResultPack.Info.MIN_NON_EDGE
                            )
                # the nested if statement here makes the combined condition mutually exclusive with the above
                if len(la2ms_r) > low_n_supporting_reads_boundary:
                    med_r = median(la2ms_r)
                    mad_r = median(map(lambda x: abs(x - med_r), la2ms_r))
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
                                | cls.ResultPack.Info.MIN_NON_EDGE
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

                    if not cond_edge_clustering:
                        info_bits |= cls.ResultPack.Info.EDGE_CLUSTERING
                    if not cond_both_strand_distrib_both:
                        info_bits |= cls.ResultPack.Info.BOTH_STRAND_DISTRIB_BOTH
                    if not cond_both_strand_distrib_one:
                        info_bits |= cls.ResultPack.Info.BOTH_STRAND_DISTRIB_ONE

                    # if no pass conditions are satisfied
                    if info_bits == (
                        cls.ResultPack.Info.EDGE_CLUSTERING
                        | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_BOTH
                        | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_ONE
                    ):
                        outcome = TestOutcomes.VARIANT_FAIL
                    if info_bits == cls.ResultPack.Info.NODATA:  # satisfied all conditions
                        info_bits |= (
                            cls.ResultPack.Info.EDGE_CLUSTERING
                            | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_BOTH
                            | cls.ResultPack.Info.BOTH_STRAND_DISTRIB_ONE
                        )

                reads_not_near_edge = (len(near_start_f) - sum(near_start_f)) + (
                    len(near_start_r) - sum(near_start_r)
                )
                if not reads_not_near_edge > min_non_edge_reads:
                    outcome = TestOutcomes.VARIANT_FAIL
                    info_bits |= cls.ResultPack.Info.MIN_NON_EDGE

                assert strand is not None  # appease type checker

                result = cls.ResultPack(outcome, strand, info_bits)

        return result


class ProportionBasedTest:
    """
    Pass/Fail a variant on proportion of supporting reads with/without a given property.

    Used by the LQF flag test to test proportion of stutter duplicate and low quality supporting reads.
    Used by the DVF flag test to test proportion of stutter duplicate reads.
    """

    @dataclass
    class ResultPack:
        class Info(Flag):
            NODATA = 0
            INSUFFCIENT_READS = 1
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
            result = cls.ResultPack(TestOutcomes.NA, cls.ResultPack.Info.INSUFFCIENT_READS, 0)
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


# END SECTION --- Variant Testing
