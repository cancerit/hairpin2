# TODO: give example in docstring?
"""
A set of Parent Abstract Base Classes that together completely describe the expected implementation of a scientific test
for recording FILTER entry for a variant record in a VCF

Heavily leverages static type checking in addition to runtime guarantees. As such, implementation of a new filter is quite verbose.
The payoff for that verbosity is:
    - complete static and runtime checking of filter validity, which is helpful when describing many filters, and complex tests
    - strong guarantees about the implementation of any given filter without needing to know about the underlying logic, making filters very easy to use once defined
"""
from abc import ABC, abstractmethod
from pydantic import BaseModel, ConfigDict
from pydantic.dataclasses import dataclass
from typing import Generic, TypeVar, Any, override
from collections.abc import Collection, Mapping
from pysam import AlignedSegment
from enum import IntEnum 
# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportUnsafeMultipleInheritance=false
# pyright: reportIncompatibleVariableOverride=false


CodeEnum_T = TypeVar("CodeEnum_T", bound=IntEnum, covariant=True)
class FilterResult(BaseModel, Generic[CodeEnum_T], ABC):
    # TODO: pydantic now obsoletes the Codes classvar, document
    """
    Parent ABC class defining the implementation that must be followed by subclasses intending to hold results of running a `FilterTester.test()` on a variant
    All filters should use subclasses that inherit from this class to hold their results.
    Defines and guarantees basic properties that must be shared by all filter results:
        - a 'getinfo' method for returning a string to report to the INFO field of the VCF record - a basic default is provided, but
          subclasses probably want to override this.
        - A basic set of instance variables, which may be extended in a subclass as necessary:
            - name, a string id for the filter to be used in the VCF FILTER field.
            - a flag, a boolean indicating if the filter is True/False, or None if untested.
            - a code, an integer code from a set of possibilities, indicating the basis on which the test has returned True/False, or None if untested.
        - A class variable, Codes, holding an IntEnum describing the set of possibilites used for the code instance variable.

    This class is generic over T, an IntEnum, such that when a subclass is made it must be made with reference to a specific IntEnum -
    this is the static equivalent of the runtime checking by Codes. When subclassing, create an appropriate enum of codes and define
    the subclass like so:
    ```
    def XYZCodes(IntEnum):
        ...
    ... # decorators if needed
    class XYZResult(FilterResult[XYZCodes]):
        Codes: ClassVar[type[XYZCodes]] = XYZCodes
        ... # the rest of the class body, e.g. further instance variables to be associated with a particular result
    ```
    """
    name: str  # for VCF FILTER field
    flag: bool | None
    code: CodeEnum_T | None

    model_config: ConfigDict = ConfigDict(frozen=True, strict=True)

    @override
    def model_post_init(self, __context: Any) -> None:
        if self.flag is not None and not self.code:
            raise ValueError('If flag is set a code must be provided')

    @abstractmethod
    def getinfo(self) -> str | None:
        """
        Return basic filter info in a string formatted for use in the VCF INFO field - "<flag>|<code>".

        Each filter must return INFO as it should be formatted for the VCF INFO field, or None if not applicable.
        Subclasses must override this method to return more specific info.
        """


@dataclass(slots=True, frozen=True)
class FilterParams:
    pass
ReadCollection_T = TypeVar("ReadCollection_T", Collection[AlignedSegment], Mapping[Any, Collection[AlignedSegment]])
FilterParams_T = TypeVar("FilterParams_T", bound=FilterParams)
FilterResult_T = TypeVar("FilterResult_T", bound=FilterResult[IntEnum], covariant=True)  # covariant such that a test method that returns a subtype of FilterResult[IntEnum] is accepted where FilterResult[IntEnum] (or FilterResult_T) is expected


class FilterTester(BaseModel, Generic[ReadCollection_T, FilterParams_T, FilterResult_T], ABC):
    """
    Parent ABC class to be inherited from when implementing a filter test on read data for a variant.
    Contains a single abstract class method, `test()`, that must be overridden by subclasses for inidvidual filters.
    """
    fixed_params: FilterParams_T

    model_config: ConfigDict = ConfigDict(frozen=True, strict=True)

    # TODO: add explanation as to how the class is generic to the docstring
    @abstractmethod
    def test(
        self,
        *args: Any,
        **kwargs: Any
    ) -> tuple[ReadCollection_T, FilterResult_T]:
        """
        Each filter must define a test method via override, respecting the method signature.
        """
