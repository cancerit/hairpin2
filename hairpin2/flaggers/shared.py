from pydantic.dataclasses import dataclass
from hairpin2.abstractflaggers import VarParams

@dataclass(slots=True)
class AltVarParams(VarParams):
    alt: str
