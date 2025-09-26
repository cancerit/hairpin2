from array import array
from dataclasses import dataclass
from enum import Flag
from statistics import mean, median, stdev
from typing import Any, cast
from collections.abc import Iterable, Sequence
from htsflow.structures import ExtendedRead, TestOutcomes, mark_read, read_has_mark, record_operation

from hairpin2.const import MutTypes, Strand, Tags, ValidatorFlags
from hairpin2.utils.ref2seq import ref2querypos, ref_end_via_cigar


# Used to tag reads as supporting a variant
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis
def check_read_supporting(
    read: ExtendedRead,
    mut_type: MutTypes,
    alt: str,
    vcf_start: int,
    vcf_stop: int,
    mark: bool = True
) -> bool:
    """
    pc8's Support Assessor

    Confirm whether a read supports a variant, and optionally mark supporting reads
    """
    if mut_type not in MutTypes:
        raise ValueError(
            f'Unsupported mutation type {mut_type} provided to mut_type argument - supports {MutTypes} only'
        )

    support_flag = ValidatorFlags.CLEAR
    try:
        mate_cig = str(read.get_tag('MC'))
    except KeyError:
        mate_cig = None
    else:
        if (not mate_cig[0].isdigit() or
            not all([(c.isalnum() or c == '=') for c in mate_cig]) or
                len(mate_cig) < 2):
            mate_cig = None
    if any(flg is None for flg in
            [read.reference_end,
                read.query_sequence,
                read.query_qualities,
                read.query_alignment_qualities,
                read.cigarstring,
                read.cigartuples,
                mate_cig]):
        support_flag |= ValidatorFlags.READ_FIELDS_MISSING
    elif mut_type in [MutTypes.SUB, MutTypes.INS]:
        try:
            mut_pos = ref2querypos(read, vcf_start)
        except ValueError:
            support_flag |= ValidatorFlags.NOT_ALIGNED
        else:
            if mut_type == MutTypes.SUB:
                if cast(str, read.query_sequence)[mut_pos:mut_pos + len(alt)] != alt:
                    support_flag |= ValidatorFlags.NOT_ALT
            if mut_type == MutTypes.INS:  # INS - mut_pos is position immediately before insertion
                if mut_pos + len(alt) > read.query_length:
                    support_flag |= ValidatorFlags.SHORT
                else:
                    mut_alns = [(q, r)
                                for q, r
                                in read.get_aligned_pairs()
                                if q in range(mut_pos + 1, mut_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        support_flag |= ValidatorFlags.BAD_OP
                    if cast(str, read.query_sequence)[mut_pos + 1:mut_pos + len(alt) + 1] != alt:
                        support_flag |= ValidatorFlags.NOT_ALT
    elif mut_type == MutTypes.DEL:
        rng = list(range(vcf_start, vcf_stop + 1))
        mut_alns = [q
                    for q, r
                    in read.get_aligned_pairs()
                    if r in rng]
        if len(mut_alns) != len(rng):
            support_flag |= ValidatorFlags.SHORT
        if (any([x is not None for x in mut_alns[1:-1]]) or
            any([x is None for x in [mut_alns[0], mut_alns[-1]]])):
                support_flag |= ValidatorFlags.BAD_OP

    # ValidatorFlag not returned, but kept in case useful in future
    if support_flag == ValidatorFlags.CLEAR:
        support = True
        if mark:
            mark_read(read, Tags.SUPPORT_TAG)
    else:
        support = False

    record_operation(read, 'mark-support')
    return support


# The result of this process is used to mark reads as low quality.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis.
# The fraction of reads tagged as low quality, or as duplicates
# by the duplicate marking function below, is used to call the
# LQF flag on variants
def check_low_qual_read(
    read: ExtendedRead,
    vcf_start: int,
    alt: str,
    mut_type: MutTypes,
    min_basequal: int,
    min_mapqual: int,
    min_clipqual: int,
    mark: bool = True
) -> bool:
    """
    pc8's Quality Assessor

    Determine if a read is of low quality, and optionally mark low quality reads
    """
    qual_flag = ValidatorFlags.CLEAR  # 0 - evaluates false


    if not (read.flag & 0x2) or read.flag & 0xE00:
        qual_flag |= ValidatorFlags.FLAG

    if read.mapping_quality < min_mapqual:
        qual_flag |= ValidatorFlags.MAPQUAL

    if ('S' in read.cigarstring and  # pyright: ignore[reportOperatorIssue]
            mean(read.query_alignment_qualities) < min_clipqual):  # pyright: ignore[reportUnknownMemberType, reportArgumentType]
        qual_flag |= ValidatorFlags.CLIPQUAL

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
                qual_flag |= ValidatorFlags.BASEQUAL

    if qual_flag == ValidatorFlags.CLEAR:
        low_qual = False
    else:
        low_qual = True
        if mark:
            mark_read(read, Tags.LOW_QUAL_TAG)

    record_operation(read, 'mark-low-qual')
    return low_qual


# Used to tag reads as an overlapping mate.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis
def check_fragment_overlap(
    read: ExtendedRead,
    vcf_start: int,
    mark: bool = True
):
    """
    pc8's Overlap Assessor

    Determine if a read is an overlapping mate for this position, and optionally mark overlapping mates
    """
    # avoid analysing both read1 and mate if they both cover the variant
    overlap = False

    mate_cig = str(read.get_tag('MC'))  # will error if no tag

    # NOTE: introduces strand bias!!
    if read.flag & 0x80:  # if second in pair
        read_range = range(
            read.reference_start,
            read.reference_end  # pyright: ignore[reportArgumentType]
        )
        mate_range = range(
            read.next_reference_start,
            ref_end_via_cigar(
                mate_cig,
                read.next_reference_start
            )
        )
        overlapping_positions = set(read_range).intersection(mate_range)
        if vcf_start in overlapping_positions:
            overlap = True

    if overlap and mark:
        mark_read(read, Tags.OVERLAP_TAG)

    record_operation(read, 'mark-overlap')

    return overlap


# The result of this process is used to mark reads as duplicates.
# This tag can then be used by other processes
# to include/exclude reads with that tag from their analysis.
# The fraction of reads tagged as duplicates by this process is
# used to assess a variant for the DVF flag
# and, partly (with the low quality tag), the LQF flag
def tag_reads_as_duplicates(
    reads: Iterable[ExtendedRead],
    duplication_window_size: int = 6
):
    """
    pc8's Stutter Duplicates Assessor

    Tag reads as duplicates which have shifted endpoints due to PCR stutter
    """
    if not reads:
        return

    # prep data
    sample_pair_endpoints: list[tuple[ExtendedRead, tuple[int, int, int, int]]] = []
    for read in reads:
        record_operation(read, 'mark-duplicates')
        next_ref_end = ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
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
        # record_operation(read, 'mark-duplicates')
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
        if all([x <= duplication_window_size for x in max_diff_per_comparison]):
            # then the read pair being examined is a duplicate of the others in the pool
            dup_endpoint_comparison_pool.append(testing_endpoints)  # add it to the comparison pool
            mark_read(read, Tags.STUTTER_DUP_TAG)
        else:
            # read at i is not dup of reads in dup_endpoint_test_pool
            # start again, test read at i against reads subsequent from i in ends_sorted
            dup_endpoint_comparison_pool = [testing_endpoints]



# Alignment Score Test -- ALF flag
# flags must be powers of 2
class ASConds(Flag):
    NODATA = 0
    INSUFFICIENT_READS = 1
    INSUFFICIENT_AS_TAGS = 2
    ON_THRESHOLD = 4

@dataclass
class ASResultPack:
    outcome: TestOutcomes
    avg_as: float | None
    reason: ASConds

def test_alignment_score(
    reads: Iterable[ExtendedRead],
    avg_AS_threshold: float
):
    """
    pc8's Alignment Score Assessor

    Pass/Fail a variant on the average alignment score of supporting reads
    """
    if not reads:
        result = ASResultPack(
            TestOutcomes.NA,
            None,
            ASConds.INSUFFICIENT_READS
        )
    else:
        aln_scores: list[float] = []

        for read in reads:
            try:
                aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
            except KeyError:
                pass
        if len(aln_scores) != 0:
            avg_as = median(aln_scores)
            if avg_as <= avg_AS_threshold:
                outcome = TestOutcomes.VARIANT_FAIL
            else:
                outcome = TestOutcomes.VARIANT_PASS

            result = ASResultPack(
                outcome,
                avg_as,
                ASConds.ON_THRESHOLD
            )
        else:
            result = ASResultPack(
                TestOutcomes.NA,
                None,
                ASConds.INSUFFICIENT_AS_TAGS
            )

    return result



# Anomalous Distribution Test -- ADF flag
# flags must be powers of 2
class ADConds(Flag):
    NODATA = 0
    NO_TESTABLE_READS = 1
    INSUFFICIENT_READS = 2
    EDGE_CLUSTERING = 4
    ONE_STRAND_DISTRIB = 8
    BOTH_STRAND_DISTRIB_BOTH = 16
    BOTH_STRAND_DISTRIB_ONE = 32
    MIN_NON_EDGE = 64


@dataclass
class ADResultPack:
    outcome: TestOutcomes
    strand: Strand
    reason: ADConds

# NOTE:
# per paper, can set hairpin for mutations distant from alignment start
# in the case where both strands have sufficient supporting reads
# discuss with Peter
def test_anomalous_distribution(
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
    min_non_edge_reads: int
) -> ADResultPack:
    """
    Pass/Fail a variant on the distribution of the variant position on the supporting reads
    
    Uses the conditions given by Ellis et al. in doi.org/10.1038/s41596-020-00437-6, and one additional condition devised by al35

    The full text of the conditional is reproduced below. There is a point of ambiguity in the original conditional.
    The interpretation that this tool has opted for is indicated by [], and is expanded upon subsequently.

    ---
    For each variant, if the number of variant-supporting reads determined in Step 59 is low (i.e., 0–1
    reads) for one strand, follow Option A. For each variant, if both strands have sufficient variant-
    supported reads (i.e., ≥2 reads), follow Option B.
    (A) Low number of variant-supporting reads on one strand
        (i) For each variant, if one strand had too few variant-supporting reads, the other strand must
        conform to:
            ● Fewer than 90% of variant-supporting reads [ON THE STRAND] have the variant located within the first 15%
            of the read measured from the alignment start position.
            ● MAD >0 and s.d. >4 for that strand.
    (B) Sufficient variant-supporting reads on both strands
        (i) For each variant, if both strands have sufficient variant-supporting reads (i.e., ≥2 reads),
        then one of the following must be true:
            ● Fewer than 90% of variant-supporting reads should have the variant located within the
            first 15% of the read measured from the alignment start position.
            ● MAD >2 and s.d. >2 for both strands.
            ● MAD >1 and s.d. >10 for one strand (i.e., strong evidence of high variability in variant
            position in variant-supporting reads).
    ---

    The point of ambiguity is whether, on path A, to include the single read from the strand which does not sufficiently
    support the variant in the test of positional distribution across the supporting reads.
    The present interpretation is that the test should be applied only to the reads on the "other strand",
    since the phrasing "the other strand must conform to" implies the exclusion of the single read on the low-support
    strand.

    The additional condition to be checked is simply whether at least N supporting reads express away from the read
    edge. This condition is experimental at this time and can be disabled by setting min_non_edge_reads to 0
    """

    
    # NOTE/TODO: test across all samples has always been done,
    # but is it really what we want to do?
    # Should it be a user choice whether to execute per sample?

    if not reads:
        result = ADResultPack(
            TestOutcomes.NA,
            Strand.BOTH,
            ADConds.NO_TESTABLE_READS,
        )
        # fresult = ResultADF(
        #     variant_flagged=TO.NA,
        #     info_flag=InfoFlagsADF.NO_TESTABLE_READS,
        #     alt=run_params.alt,
        #     reads_seen=0,
        #     strand='BOTH'
        # )
    else:
        # *l*engths of *a*lignment starts *to* *m*utant query positions
        la2ms_f: list[int] = []
        la2ms_r: list[int] = []
        near_start_f: list[bool] = []
        near_start_r: list[bool] = []

        for read in reads:
            try:
                mut_qpos = ref2querypos(read, record_start)  # TODO: this is additive processing and should be eventually separated as such
            except ValueError:
                raise ValueError(f'read {read.query_name} does not cover variant')

            if read.flag & 0x10:
                # +1 to include last base in length
                la2m = read.query_alignment_end - mut_qpos + 1
                near_start_r.append(
                    (la2m / read.query_alignment_length) <= edge_definition
                )
                la2ms_r.append(la2m)
            else:
                la2m = mut_qpos - read.query_alignment_start + 1
                near_start_f.append(
                    (la2m / read.query_alignment_length) <= edge_definition
                )
                la2ms_f.append(la2m)

        # hairpin conditions from Ellis et al. 2020, Nature Protocols
        # sometimes reported as 2021
        if len(la2ms_f) <= low_n_supporting_reads_boundary and len(la2ms_r) <= low_n_supporting_reads_boundary:
            result = ADResultPack(
                TestOutcomes.NA,
                Strand.BOTH,
                ADConds.INSUFFICIENT_READS
            )
        else:
            already_passed = ADConds.NO_TESTABLE_READS | ADConds.INSUFFICIENT_READS
            strand = None
            outcome =  TestOutcomes.VARIANT_PASS
            info_bits = ADConds.NODATA
            if len(la2ms_f) > low_n_supporting_reads_boundary:  # if true, calculate stats
                med_f = median(la2ms_f)  # range calculation replaced with true MAD calc (for r strand also)
                mad_f = median(map(lambda x: abs(x - med_f), la2ms_f))
                sd_f = stdev(la2ms_f)
                if len(la2ms_r) <= low_n_supporting_reads_boundary:  # if also this, test on this path
                    strand = Strand.F
                    cond_edge_clustering = (sum(near_start_f) / len(near_start_f)) < edge_clustering_threshold
                    cond_one_strand_mad = mad_f > min_MAD_one_strand
                    cond_one_strand_sd = sd_f > min_sd_one_strand

                    # this could be more concise, but I think this is the most readable approach I've found
                    if not cond_edge_clustering:
                        info_bits |= ADConds.EDGE_CLUSTERING
                        outcome = TestOutcomes.VARIANT_FAIL
                    if not (cond_one_strand_mad and cond_one_strand_sd):
                        info_bits |= ADConds.ONE_STRAND_DISTRIB
                        outcome = TestOutcomes.VARIANT_FAIL

                    if outcome == TestOutcomes.VARIANT_PASS:
                        info_bits = ADConds.ONE_STRAND_DISTRIB | ADConds.EDGE_CLUSTERING | ADConds.MIN_NON_EDGE | already_passed
            # the nested if statement here makes the combined condition mutually exclusive with the above
            if len(la2ms_r) > low_n_supporting_reads_boundary:
                med_r = median(la2ms_r)
                mad_r = median(map(lambda x: abs(x - med_r), la2ms_r))
                sd_r = stdev(la2ms_r)
                if len(la2ms_f) <= low_n_supporting_reads_boundary:
                    strand = Strand.R
                    outcome = TestOutcomes.VARIANT_PASS
                    info_bits = ADConds.NODATA
                    cond_edge_clustering = (sum(near_start_r) / len(near_start_r)) < edge_clustering_threshold
                    cond_one_strand_mad = mad_r > min_MAD_one_strand
                    cond_one_strand_sd = sd_r > min_sd_one_strand

                    if not cond_edge_clustering:
                        info_bits |= ADConds.EDGE_CLUSTERING
                        outcome = TestOutcomes.VARIANT_FAIL
                    if not (cond_one_strand_mad and cond_one_strand_sd):
                        info_bits |= ADConds.ONE_STRAND_DISTRIB
                        outcome = TestOutcomes.VARIANT_FAIL

                    if outcome == TestOutcomes.VARIANT_PASS:
                        info_bits = ADConds.ONE_STRAND_DISTRIB | ADConds.EDGE_CLUSTERING | ADConds.MIN_NON_EDGE | already_passed
            if len(la2ms_f) > low_n_supporting_reads_boundary and len(la2ms_r) > low_n_supporting_reads_boundary:
                strand = Strand.BOTH
                outcome = TestOutcomes.VARIANT_PASS
                info_bits = ADConds.NODATA
                frac_lt_thresh = (
                    sum(near_start_f + near_start_r)
                    / (len(near_start_f) + len(near_start_r))
                )

                cond_edge_clustering = frac_lt_thresh < edge_clustering_threshold
                cond_both_strand_distrib_both = (
                    mad_f > min_MAD_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    mad_r > min_MAD_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    sd_f > min_sd_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    sd_r > min_sd_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                )
                cond_both_strand_distrib_one = (
                    (mad_f > min_MAD_both_strand_strong and sd_f > min_sd_both_strand_strong) or  # pyright: ignore[reportPossiblyUnboundVariable]
                    (mad_r > min_MAD_both_strand_strong and sd_r > min_sd_both_strand_strong)  # pyright: ignore[reportPossiblyUnboundVariable]
                )

                if not cond_edge_clustering:
                    info_bits |= ADConds.EDGE_CLUSTERING
                if not cond_both_strand_distrib_both:
                    info_bits |= ADConds.BOTH_STRAND_DISTRIB_BOTH
                if not cond_both_strand_distrib_one:
                    info_bits |= ADConds.BOTH_STRAND_DISTRIB_ONE

                # if no pass conditions are satisfied
                if info_bits == (ADConds.EDGE_CLUSTERING | ADConds.BOTH_STRAND_DISTRIB_BOTH | ADConds.BOTH_STRAND_DISTRIB_ONE):
                    outcome = TestOutcomes.VARIANT_FAIL
                if info_bits == ADConds.NODATA:  # satisfied all conditions
                    info_bits |= ADConds.EDGE_CLUSTERING | ADConds.BOTH_STRAND_DISTRIB_BOTH | ADConds.BOTH_STRAND_DISTRIB_ONE | already_passed

            reads_not_near_edge = (len(near_start_f) - sum(near_start_f)) + (len(near_start_r) - sum(near_start_r))
            if not reads_not_near_edge > min_non_edge_reads:
                outcome = TestOutcomes.VARIANT_FAIL
                info_bits |= ADConds.MIN_NON_EDGE

            assert strand is not None  # appease type checker

            result = ADResultPack(
                outcome,
                strand,
                info_bits
            )

    return result



# --- ASSESS PROPORTION OF READS WITH A GIVEN PROPERTY ---
# used by LQF (checks duplicates, low qual)
# and DVF (checks duplicates)
class PropConds(Flag):
    NODATA = 0
    INSUFFCIENT_READS = 1
    THRESHOLD = 2
    MIN_PASS = 4

@dataclass
class PropResultPack:
    outcome: TestOutcomes
    reason: PropConds
    prop_loss: float

def test_proportion_with_tag(
    reads: Sequence[ExtendedRead],
    tags_to_check: Iterable[Tags],
    read_loss_threshold: float,
    min_without: int = 0
) -> PropResultPack:
    if not reads:
        result = PropResultPack(
            TestOutcomes.NA,
            PropConds.INSUFFCIENT_READS,
            0
        )
    else:
        info_bits = PropConds.NODATA
        outcome = TestOutcomes.VARIANT_PASS
        ntotal = len(reads)

        nwith = 0
        for read in reads:
            nwith += int(bool(sum(read_has_mark(read, tag) for tag in tags_to_check)))
        nwithout = ntotal - nwith
        loss_ratio= nwith / ntotal  # you could alternately formulate this as nwith / nwithout...
        if loss_ratio > read_loss_threshold:
            info_bits |= PropConds.THRESHOLD
            outcome = TestOutcomes.VARIANT_FAIL
        if nwithout < min_without:
            info_bits |= PropConds.MIN_PASS
            outcome = TestOutcomes.VARIANT_FAIL

        if outcome == TestOutcomes.VARIANT_PASS:
            info_bits = ~PropConds(0)  # satisfied all conditions

        result = PropResultPack(
            outcome,
            info_bits,
            loss_ratio
        )

    return result
