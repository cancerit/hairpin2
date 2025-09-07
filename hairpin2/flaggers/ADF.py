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
import hairpin2.abstractflaggers as haf
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2 import ref2seq as r2s
from typing import override
from enum import IntEnum, auto
from statistics import median, stdev

FLAG_NAME = "ADF"


class CodesADF(IntEnum):
    INSUFFICIENT_READS = 0
    SIXTYAI = auto()
    SIXTYBI = auto()


class ResultADF(
    haf.FlagResult,
    flag_name=FLAG_NAME,
    result_codes=tuple(CodesADF)
):
    alt: str

    @override
    def getinfo(self) -> str:
        return f'{self.alt}|{self.flag}|{self.code}'


class FixedParamsADF(haf.FixedParams):
    edge_definition: float = 0.15  # relative proportion, by percentage, of a read to be considered 'the edge'
    edge_clustering_threshold: float = 0.9  # percentage threshold
    min_MAD_one_strand: int = 0  # exclusive (and subsequent params)
    min_sd_one_strand: float = 4
    min_MAD_both_strand_weak: int = 2
    min_sd_both_strand_weak: float = 2
    min_MAD_both_strand_strong: int = 1
    min_sd_both_strand_strong: float = 10
    min_reads: int = 1  # inclusive


# def prefilter(
#     self
# ):
#     filtered_reads: dict[Any, list[AlignedSegment]] = {}
#     for sample_key, reads in self.run_params.reads.items():
#         passed_reads: list[AlignedSegment] = []
#         for read in reads:
#             invalid = qc_read(
#                 read,
#                 self.run_params.record.start,
#                 self.run_params.record.stop,
#                 self.run_params.alt,
#                 self.run_params.mut_type,
#                 self.prefilter_params.min_baseq,
#                 self.prefilter_params.min_mapq,
#                 self.prefilter_params.min_avg_clipq,
                
#             )
#             if not invalid:
#                 passed_reads.append(read)
#         filtered_reads[sample_key] = passed_reads
#     return filtered_reads


# NOTE/TODO:
# per paper, can set hairpin for mutations distant from alignment start
# in the case where both strands have sufficient supporting reads
# discuss with Peter
def test_anomalous_distribution(
    run_params: RunParamsShared,
    fixed_params: FixedParamsADF
) -> ResultADF:
    # NOTE/TODO: test across all samples has always been done,
    # but is it really what we want to do?
    # Should it be a user choice whether to execute per sample?

    if len(run_params.reads.all) < 1:
        fresult = ResultADF(
            flag=None,
            code=CodesADF.INSUFFICIENT_READS,
            alt=run_params.alt
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
        if len(la2ms_f) <= fixed_params.min_reads and len(la2ms_r) <= fixed_params.min_reads:
            fresult = ResultADF(
                flag=None,
                code=CodesADF.INSUFFICIENT_READS,
                alt=run_params.alt
            )
        else:
            if len(la2ms_f) > fixed_params.min_reads:  # if true, calculate stats
                med_f = median(la2ms_f)  # range calculation replaced with true MAD calc (for r strand also)
                mad_f = median(map(lambda x: abs(x - med_f), la2ms_f))
                sd_f = stdev(la2ms_f)
                if len(la2ms_r) <= fixed_params.min_reads:  # if also this, test
                    if (((sum(near_start_f) / len(near_start_f)) < fixed_params.edge_clustering_threshold) and
                        mad_f > fixed_params.min_MAD_one_strand and
                            sd_f > fixed_params.min_sd_one_strand):
                        code = CodesADF.SIXTYAI  # 60A(i)
                        flag = False
                    else:
                        code = CodesADF.SIXTYAI
                        flag = True
            # the nested if statement here makes the combined condition mutually exclusive with the above
            if len(la2ms_r) > fixed_params.min_reads:
                med_r = median(la2ms_r)
                mad_r = median(map(lambda x: abs(x - med_r), la2ms_r))
                sd_r = stdev(la2ms_r)
                if len(la2ms_f) <= fixed_params.min_reads:
                    if (((sum(near_start_r) / len(near_start_r)) < fixed_params.edge_clustering_threshold) and
                        mad_r > fixed_params.min_MAD_one_strand and
                            sd_r > fixed_params.min_sd_one_strand):
                        code = CodesADF.SIXTYAI
                        flag = False
                    else:
                        code = CodesADF.SIXTYAI
                        flag = True
            if len(la2ms_f) > fixed_params.min_reads and len(la2ms_r) > fixed_params.min_reads:
                frac_lt_thresh = (sum(near_start_f + near_start_r)
                                  / (len(near_start_f) + len(near_start_r)))
                if (frac_lt_thresh < fixed_params.edge_clustering_threshold or
                    (mad_f > fixed_params.min_MAD_both_strand_weak and mad_r > fixed_params.min_MAD_both_strand_weak and sd_f > fixed_params.min_sd_both_strand_weak and sd_r > fixed_params.min_sd_both_strand_weak) or  # pyright: ignore[reportPossiblyUnboundVariable]
                    (mad_f > fixed_params.min_MAD_both_strand_strong and sd_f > fixed_params.min_sd_both_strand_strong) or  # pyright: ignore[reportPossiblyUnboundVariable]
                        (mad_r > fixed_params.min_MAD_both_strand_strong and sd_r > fixed_params.min_sd_both_strand_strong)):  # pyright: ignore[reportPossiblyUnboundVariable]
                    code = CodesADF.SIXTYBI  # 60B(i)
                    flag = False
                else:
                    code = CodesADF.SIXTYBI
                    flag = True
            fresult = ResultADF(
                flag=flag,  # pyright: ignore[reportPossiblyUnboundVariable]
                code=code,  # pyright: ignore[reportPossiblyUnboundVariable]
                alt=run_params.alt
            )

    return fresult


# require support, exclude stutter dups, overlapping second pair member, low quality
@haf.require_read_properties(require_tags=['zS'], exclude_tags=['zD', 'zO', 'zQ'])  # TODO guard against accidentally providing just a string - which is an iterable[str]!!
@haf.variant_flagger(flag_name=FLAG_NAME, flagger_param_class=FixedParamsADF, flagger_func=test_anomalous_distribution, result_type=ResultADF)
class FlaggerADF(
    haf.ReadAwareProcess
):
    # TODO: docstring
    """
    Anomalous Distribution Filter based on hairpin filtering algorthim described in Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6) 
    """
