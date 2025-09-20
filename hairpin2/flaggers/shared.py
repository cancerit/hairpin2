from typing import Literal

from hairpin2.abstractions.process_params import RunParams


class RunParamsShared(RunParams):
    alt: str
    mut_type: Literal['S', 'D', 'I']
