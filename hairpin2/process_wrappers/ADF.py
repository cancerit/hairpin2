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
from typing import Annotated, override

from pydantic import BeforeValidator

from hairpin2.const import FlaggerNamespaces, Strand
from hairpin2.infrastructure.configure_funcs import make_variant_flagger
from hairpin2.infrastructure.constraints import bound
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.infrastructure.process_params import FixedParams
from hairpin2.infrastructure.structures import FlagResult
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import AnomalousDistributionTest


@dataclass(frozen=True)
class ResultADF(
    FlagResult,
    flag_name=FlaggerNamespaces.ANOMALOUS_DISTRIBUTION,
    info_enum=AnomalousDistributionTest.ResultPack.Info,
):
    alt: str
    reads_seen: int
    strand: Strand

    @override
    def getinfo(self) -> str:
        info_bits = hex(self.info_flag.value if self.info_flag is not None else 0)
        return f"{self.alt}|{self.variant_flagged.value}|{info_bits}|{self.strand.value}|{self.reads_seen}"


def validate_positive(val: int | float):
    return bound(val, 0)


class FixedParamsADF(FixedParams):
    # relative proportion, by percentage, of a read to be considered 'the edge'
    edge_definition: Annotated[float, BeforeValidator(lambda x: bound(x, 0.0, 1.0))]
    # percentage threshold
    edge_clustering_threshold: Annotated[float, BeforeValidator(lambda x: bound(x, 0.0, 1.0))]
    # exclusive (and subsequent params)
    min_MAD_one_strand: Annotated[float, BeforeValidator(validate_positive)]
    min_sd_one_strand: Annotated[float, BeforeValidator(validate_positive)]
    min_MAD_both_strand_weak: Annotated[float, BeforeValidator(validate_positive)]
    min_sd_both_strand_weak: Annotated[float, BeforeValidator(validate_positive)]
    min_MAD_both_strand_strong: Annotated[float, BeforeValidator(validate_positive)]
    min_sd_both_strand_strong: Annotated[float, BeforeValidator(validate_positive)]
    min_non_edge_reads: Annotated[int, BeforeValidator(validate_positive)]
    # inclusive
    low_n_supporting_reads_boundary: Annotated[int, BeforeValidator(validate_positive)]


def test_adf(run_params: RunParamsShared, fixed_params: FixedParamsADF):
    result = AnomalousDistributionTest.test_variant_reads(
        run_params.reads.all,
        run_params.record.start,
        fixed_params.edge_definition,
        fixed_params.edge_clustering_threshold,
        fixed_params.min_MAD_one_strand,
        fixed_params.min_sd_one_strand,
        fixed_params.min_MAD_both_strand_weak,
        fixed_params.min_sd_both_strand_weak,
        fixed_params.min_MAD_both_strand_strong,
        fixed_params.min_sd_both_strand_strong,
        fixed_params.low_n_supporting_reads_boundary,
        fixed_params.min_non_edge_reads,
    )

    flag = ResultADF(
        result.outcome, result.reason, run_params.alt, len(run_params.reads.all), result.strand
    )

    return flag


# require support, exclude stutter dups, overlapping second pair member, low quality
@make_variant_flagger(
    process_namespace=FlaggerNamespaces.ANOMALOUS_DISTRIBUTION,
    flagger_param_class=FixedParamsADF,
    flagger_func=test_adf,
    result_type=ResultADF,
)
class FlaggerADF(
    ReadAwareProcess,
):
    """
    Anomalous Distribution Filter based on the hairpin filtering algorithm described in Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6)
    """
