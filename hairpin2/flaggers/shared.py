from typing import Literal
from hairpin2.abstractions.readawareproc import FixedParams, RunParams


class RunParamsShared(RunParams):
    alt: str
    mut_type: Literal['S', 'D', 'I']


class PrefilterParamsShared(FixedParams):
    min_mapping_quality: int
    min_avg_clip_quality: int
    min_base_quality: int
