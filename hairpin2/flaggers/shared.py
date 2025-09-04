from typing import Literal
from hairpin2.abstractflaggers import PrefilterParams, RunParams


class RunParamsShared(RunParams):
    alt: str
    mut_type: Literal['S', 'D', 'I']


class PrefilterParamsShared(PrefilterParams):
    min_mapq: int
    min_avg_clipq: int
    min_baseq: int
