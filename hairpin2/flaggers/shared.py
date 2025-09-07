from typing import Literal
from hairpin2.abstractflaggers import FixedParams, RunParams


class RunParamsShared(RunParams):
    alt: str
    mut_type: Literal['S', 'D', 'I']


class PrefilterParamsShared(FixedParams):
    min_mapq: int
    min_avg_clipq: int
    min_baseq: int
