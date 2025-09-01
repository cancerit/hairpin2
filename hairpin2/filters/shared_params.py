from pydantic.dataclasses import dataclass
from hairpin2.abstractfilters import VarParams

@dataclass(slots=True)
class AltVarParams(VarParams):
    alt: str
