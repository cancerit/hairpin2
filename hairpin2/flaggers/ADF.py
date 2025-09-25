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
from hairpin2.abstractions.configure_funcs import make_variant_flagger
from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.process_params import FixedParams
from hairpin2.abstractions.structures import FlagResult
from hairpin2.const import FlaggerNamespaces, TestOutcomes as TO
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.utils import ref2seq as r2s
from typing import Literal, cast, override
from enum import IntFlag, auto
from statistics import median, stdev


# BUG: codes changed and isn't documented
# TODO: flags allow for greater information packing. Use flags, and provide decomp cli tool
class InfoFlagsADF(IntFlag):
    NO_TESTABLE_READS = 1
    INSUFFICIENT_READS = auto()
    EDGE_CLUSTERING = auto()
    ONE_STRAND_DISTRIB = auto()
    BOTH_STRAND_DISTRIB_BOTH = auto()
    BOTH_STRAND_DISTRIB_ONE = auto()
    MIN_NON_EDGE = auto()


@dataclass(frozen=True)
class ResultADF(
    FlagResult,
    flag_name=FlaggerNamespaces.ANOMALOUS_DISTRIBUTION,
    info_enum=InfoFlagsADF
):
    alt: str
    strand: Literal['F', 'R', "BOTH"]

    @override
    def getinfo(self) -> str:
        return f'{self.alt}|{self.variant_flagged}|{self.info_flag}'


class FixedParamsADF(FixedParams):
    edge_definition: float  # relative proportion, by percentage, of a read to be considered 'the edge'
    edge_clustering_threshold: float  # percentage threshold
    min_MAD_one_strand: int  # exclusive (and subsequent params)
    min_sd_one_strand: float
    min_MAD_both_strand_weak: int
    min_sd_both_strand_weak: float
    min_MAD_both_strand_strong: int
    min_sd_both_strand_strong: float
    min_non_edge_reads: int
    low_n_supporting_reads_boundary: int  # inclusive


# NOTE:
# per paper, can set hairpin for mutations distant from alignment start
# in the case where both strands have sufficient supporting reads
# discuss with Peter
# TODO: hts-flow. The idea with the function signature decomp is to be able
# to write functions that can be used with no knowledge
# of anything defined by hts-flow such that they're reusable elsewhere.
# Exactly what this looks like is tbd. One could request knowledge of ExtendedRead
# or one could enforce annotation of parameters with Annotated. Try both
def test_anomalous_distribution(
    run_params: RunParamsShared,
    fixed_params: FixedParamsADF
) -> ResultADF:
    # NOTE/TODO: test across all samples has always been done,
    # but is it really what we want to do?
    # Should it be a user choice whether to execute per sample?

    if len(run_params.reads.all) < 1:
        fresult = ResultADF(
            variant_flagged=TO.NA,
            info_flag=InfoFlagsADF.NO_TESTABLE_READS,
            alt=run_params.alt,
            strand='BOTH'
        )
    else:
        # *l*engths of *a*lignment starts *to* *m*utant query positions
        la2ms_f: list[int] = []
        la2ms_r: list[int] = []
        near_start_f: list[bool] = []
        near_start_r: list[bool] = []

        for read in run_params.reads.all:
            try:
                mut_qpos = r2s.ref2querypos(read, run_params.record.start)  # TODO: this is additive processing and should be eventually separated as such
            except ValueError:
                raise ValueError(f'read {read.query_name} does not cover variant')

            if read.flag & 0x10:
                # +1 to include last base in length
                la2m = read.query_alignment_end - mut_qpos + 1
                near_start_r.append(
                    (la2m / read.query_alignment_length) <= fixed_params.edge_definition
                )
                la2ms_r.append(la2m)
            else:
                la2m = mut_qpos - read.query_alignment_start + 1
                near_start_f.append(
                    (la2m / read.query_alignment_length) <= fixed_params.edge_definition
                )
                la2ms_f.append(la2m)

       # hairpin conditions from Ellis et al. 2020, Nature Protocols
        # sometimes reported as 2021
        if len(la2ms_f) <= fixed_params.low_n_supporting_reads_boundary and len(la2ms_r) <= fixed_params.low_n_supporting_reads_boundary:
            fresult = ResultADF(
                variant_flagged=TO.NA,
                info_flag=InfoFlagsADF.INSUFFICIENT_READS,
                alt=run_params.alt,
                strand='BOTH'
            )
        else:
            strand = None
            info_bits = None
            flag = None
            if len(la2ms_f) > fixed_params.low_n_supporting_reads_boundary:  # if true, calculate stats
                med_f = median(la2ms_f)  # range calculation replaced with true MAD calc (for r strand also)
                mad_f = median(map(lambda x: abs(x - med_f), la2ms_f))
                sd_f = stdev(la2ms_f)
                if len(la2ms_r) <= fixed_params.low_n_supporting_reads_boundary:  # if also this, test on this path
                    strand = 'F'
                    flag =  TO.VARIANT_PASS
                    info_bits = 0
                    cond_edge_clustering = (sum(near_start_f) / len(near_start_f)) < fixed_params.edge_clustering_threshold
                    cond_one_strand_mad = mad_f > fixed_params.min_MAD_one_strand
                    cond_one_strand_sd = sd_f > fixed_params.min_sd_one_strand

                    # this could be more concise, but I think this is the most readable approach I've found
                    if not cond_edge_clustering:
                        info_bits |= InfoFlagsADF.EDGE_CLUSTERING
                        flag = TO.VARIANT_FAIL
                    if not (cond_one_strand_mad and cond_one_strand_sd):
                        info_bits |= InfoFlagsADF.ONE_STRAND_DISTRIB
                        flag = TO.VARIANT_FAIL

                    if flag == TO.VARIANT_PASS:
                        info_bits = InfoFlagsADF.ONE_STRAND_DISTRIB | InfoFlagsADF.EDGE_CLUSTERING | InfoFlagsADF.MIN_NON_EDGE
            # the nested if statement here makes the combined condition mutually exclusive with the above
            if len(la2ms_r) > fixed_params.low_n_supporting_reads_boundary:
                med_r = median(la2ms_r)
                mad_r = median(map(lambda x: abs(x - med_r), la2ms_r))
                sd_r = stdev(la2ms_r)
                if len(la2ms_f) <= fixed_params.low_n_supporting_reads_boundary:
                    strand = 'R'
                    flag = TO.VARIANT_PASS
                    info_bits = 0
                    cond_edge_clustering = (sum(near_start_r) / len(near_start_r)) < fixed_params.edge_clustering_threshold
                    cond_one_strand_mad = mad_r > fixed_params.min_MAD_one_strand
                    cond_one_strand_sd = sd_r > fixed_params.min_sd_one_strand

                    if not cond_edge_clustering:
                        info_bits |= InfoFlagsADF.EDGE_CLUSTERING
                        flag = TO.VARIANT_FAIL
                    if not (cond_one_strand_mad and cond_one_strand_sd):
                        info_bits |= InfoFlagsADF.ONE_STRAND_DISTRIB
                        flag = TO.VARIANT_FAIL

                    if flag == TO.VARIANT_PASS:
                        info_bits = InfoFlagsADF.ONE_STRAND_DISTRIB | InfoFlagsADF.EDGE_CLUSTERING | InfoFlagsADF.MIN_NON_EDGE
            if len(la2ms_f) > fixed_params.low_n_supporting_reads_boundary and len(la2ms_r) > fixed_params.low_n_supporting_reads_boundary:
                strand = 'BOTH'
                flag = TO.VARIANT_PASS
                info_bits = 0
                frac_lt_thresh = (
                    sum(near_start_f + near_start_r)
                    / (len(near_start_f) + len(near_start_r))
                )

                cond_edge_clustering = frac_lt_thresh < fixed_params.edge_clustering_threshold
                cond_both_strand_distrib_both = (
                    mad_f > fixed_params.min_MAD_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    mad_r > fixed_params.min_MAD_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    sd_f > fixed_params.min_sd_both_strand_weak and  # pyright: ignore[reportPossiblyUnboundVariable]
                    sd_r > fixed_params.min_sd_both_strand_weak  # pyright: ignore[reportPossiblyUnboundVariable]
                )
                cond_both_strand_distrib_one = (
                    (mad_f > fixed_params.min_MAD_both_strand_strong and sd_f > fixed_params.min_sd_both_strand_strong) or  # pyright: ignore[reportPossiblyUnboundVariable]
                    (mad_r > fixed_params.min_MAD_both_strand_strong and sd_r > fixed_params.min_sd_both_strand_strong)  # pyright: ignore[reportPossiblyUnboundVariable]
                )

                if not cond_edge_clustering:
                    info_bits |= InfoFlagsADF.EDGE_CLUSTERING
                if not cond_both_strand_distrib_both:
                    info_bits |= InfoFlagsADF.BOTH_STRAND_DISTRIB_BOTH
                if not cond_both_strand_distrib_one:
                    info_bits |= InfoFlagsADF.BOTH_STRAND_DISTRIB_ONE

                # if no pass conditions are satisfied
                if info_bits == (InfoFlagsADF.EDGE_CLUSTERING | InfoFlagsADF.BOTH_STRAND_DISTRIB_BOTH | InfoFlagsADF.BOTH_STRAND_DISTRIB_ONE):
                    flag = TO.VARIANT_FAIL
                if info_bits == 0:  # satisfied all conditions
                    info_bits = (InfoFlagsADF.EDGE_CLUSTERING | InfoFlagsADF.BOTH_STRAND_DISTRIB_BOTH | InfoFlagsADF.BOTH_STRAND_DISTRIB_ONE)

            # appease type checker (and guard)
            assert strand is not None
            assert info_bits is not None
            assert flag is not None

            reads_not_near_edge = (len(near_start_f) - sum(near_start_f)) + (len(near_start_r) - sum(near_start_r))
            if not reads_not_near_edge > fixed_params.min_non_edge_reads:
                flag = TO.VARIANT_FAIL
                info_bits |= InfoFlagsADF.MIN_NON_EDGE

            fresult = ResultADF(
                variant_flagged=flag,
                info_flag=cast(InfoFlagsADF, info_bits),
                alt=run_params.alt,
                strand=strand
            )

    return fresult


# require support, exclude stutter dups, overlapping second pair member, low quality
@make_variant_flagger(
    process_namespace=FlaggerNamespaces.ANOMALOUS_DISTRIBUTION,
    flagger_param_class=FixedParamsADF,
    flagger_func=test_anomalous_distribution,
    result_type=ResultADF
)
class FlaggerADF(
    ReadAwareProcess,
):
    """
    Anomalous Distribution Filter based on hairpin filtering algorthim described in Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6) 
    """
