# SECTION: Params --------------------------------
# Users define the params required for a process
# using the nice pydantic dataclass-like syntax
# and then they can use those in their process functions
# with dot syntax. nice. Also free pydantic validation

from pydantic import BaseModel, ConfigDict
from pysam import VariantRecord
from hairpin2.infrastructure.structures import ExtendedRead, ReadView


class _Params(BaseModel):
    model_config: ConfigDict = ConfigDict(  # pyright: ignore[reportIncompatibleVariableOverride]
        strict=True, frozen=True, arbitrary_types_allowed=True, extra="forbid"
    )

    """
    parent dataclass to be be inherited from to store specific fixed parameters for a particular subclass of FilterTester,
    or in other words for a particular filtering test. Using subclasses of this class for the fixed parameters provides
    type-safety and a consistent interface for implementing filters
    """
    pass


# inherit and set defaults in subclasses
class FixedParams(_Params):
    pass


# some mandatory fields
class RunParams(_Params):
    record: VariantRecord
    reads: ReadView[ExtendedRead]
