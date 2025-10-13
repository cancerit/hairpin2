from hairpin2.const import MutTypes
from hairpin2.infrastructure.process_params import RunParams


class RunParamsShared(RunParams):
    alt: str
    mut_type: MutTypes
