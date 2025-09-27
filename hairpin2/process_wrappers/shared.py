from hairpin2.infrastructure.process_params import RunParams
from hairpin2.const import MutTypes


class RunParamsShared(RunParams):
    alt: str
    mut_type: MutTypes
