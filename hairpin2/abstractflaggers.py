# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportUnsafeMultipleInheritance=false
# pyright: reportIncompatibleVariableOverride=false
"""
A set of Parent Abstract Base Classes that together completely describe the expected implementation of a scientific test
for recording FILTER entry for a variant record in a VCF

Heavily leverages static type checking in addition to runtime guarantees. As such, implementation of a new filter is quite verbose.
The payoff for that verbosity is:
    - complete static and runtime checking of filter validity, which is helpful when describing many filters, and complex tests
    - strong guarantees about the implementation of any given filter without needing to know about the underlying logic, making filters very easy to use once defined
"""
from abc import ABC, abstractmethod
from pydantic import BaseModel, ConfigDict, PrivateAttr
from pydantic.dataclasses import dataclass
from typing import Generic, TypeVar, Any, override, ClassVar
from pysam import VariantRecord
from enum import IntEnum
from hairpin2.structures import ReadView


### READ FilterTester First and work backwards! ###

CodeEnum_T = TypeVar("CodeEnum_T", bound=IntEnum, covariant=True)
class FilterResult(BaseModel, Generic[CodeEnum_T], ABC):
    """
    Parent ABC/pydantic BaseModel defining the implementation that must be followed by result subclasses - subclasses to hold results of running
    the `test()` method of a specific subclass of FilterTester on a variant (e.g. for `FooFilter.test()`, the return value should include an instance of `FooResult`).
    All filters should use subclasses that inherit from this class to hold their results, as enforced by the FilterTester ABC.
    Defines and guarantees basic properties that must be shared by all filter results:
        - a 'getinfo' method for returning a string to report to the INFO field of the VCF record - a basic default is provided, but
          subclasses probably want to override this.
        - A basic set of instance variables, which may be extended in a subclass as necessary:
            - name, a string id for the filter to be used in the VCF FILTER field.
            - a flag, a boolean indicating if the filter is True/False, or None if untested.
            - a code, an integer code from a set of possibilities, indicating the basis on which the test has returned True/False, or None if untested.

    This class is generic over T, an IntEnum, such that when a subclass is made it must be made with reference to a specific IntEnum -
    allowing for static analysis to ensure only the correct codes are set in the result instance. Since this parent class is a pydantic BaseModel,
    subclasses also have runtime validation of parameters, including these codes.

    When subclassing, create an appropriate enum of codes and define the subclass like so:
    ```
    def XYZCodes(IntEnum):
        ...
    ... # decorators if needed
    class XYZResult(FilterResult[XYZCodes]):
        ... # the rest of the class body, e.g. further instance variables to be associated with a particular result
    ```
    """
    Name: ClassVar[str]  # for VCF FILTER field
    flag: bool | None
    code: CodeEnum_T | None

    model_config: ConfigDict = ConfigDict(frozen=True, strict=True)

    @override
    def model_post_init(self, __context: Any) -> None:
        if self.flag is not None and not self.code:
            raise ValueError('If flag is set a code must be provided')

    @abstractmethod
    def getinfo(self) -> str:
        """
        Return basic filter info in a string formatted for use in the VCF INFO field - "<flag>|<code>".

        Each filter must return INFO as it should be formatted for the VCF INFO field, or None if not applicable.
        Subclasses must override this method to return more specific info.
        """


class Params:
    '''
    parent dataclass to be be inherited from to store specific fixed parameters for a particular subclass of FilterTester,
    or in other words for a particular filtering test. Using subclasses of this class for the fixed parameters provides
    type-safety and a consistent interface for implementing filters
    '''
    pass


@dataclass(frozen=True)
class FixedParams(Params):
    pass


@dataclass(
    config=ConfigDict(
        frozen=True,
        strict=True,
        arbitrary_types_allowed=True
    )
)
class VarParams(Params):
    record: VariantRecord  # makes sense right?
    reads: ReadView


FixedParams_T = TypeVar("FixedParams_T", bound=FixedParams)
VarParams_T = TypeVar("VarParams_T", bound=VarParams)
FilterResult_T = TypeVar("FilterResult_T", bound=FilterResult[IntEnum], covariant=True)  # covariant such that a test method that returns a subtype of FilterResult[IntEnum] is accepted where FilterResult[IntEnum] (or FilterResult_T) is expected


class Flagger(BaseModel, Generic[FixedParams_T, VarParams_T, FilterResult_T], ABC):
    # TODO: update docstring
    """
    Parent ABC/pydantic BaseModel to be inherited from when implementing a filter test on read data for a variant.
    Contains a single abstract class method, `test()`, that must be overridden by subclasses for inidvidual filters.

    The class is Generic over type variables ReadCollection_T (Collection[AlignedSegment]/Mapping[Any, Collection[AlignedSegment]]),
    FilterParams_T, and FilterResult_t. This means that on implementing a concrete filter by subclassing this class,
    you must select the specific types to which the filter pertains. For example, if implementing filter FooFilter,
    the defintion might look something like:
    `class FooFilter(FilterTester[list[AlignedSegment], FooParams, FooResult])`

    As a result:
        - static analysis will insist that the `test()` method of FooFilter returns tuple[list[AlignedSegment], FooResult]
        - static analysis will insist instances of FooFilter pass the correct params upon instantiation
        - since FilterTester is a pydantic BaseModel,
          FooFilter will also validate that the correct FilterParams have been passed at runtime

    In total, using this parent class ensures concrete filters obey strict rules, and therefore:
        - filters are more difficult to implement incorrectly
        - filters have a sufficiently consistent interface such that adding extra filters into main becomes trivial
    """
    fixed_params: FixedParams_T
    _var_params: VarParams_T | None = PrivateAttr(None)
    _executed: bool = PrivateAttr(False)

    model_config: ConfigDict = ConfigDict(
        frozen=True,
        strict=True,
    )

    # read only
    @property
    def var_params(self) -> VarParams_T:
        if self._var_params is None:
            raise AttributeError("var_params not set, cannot access")
        return self._var_params

    def prime(
        self,
        var_params: VarParams_T,
        overwrite: bool = False
    ):
        if self._executed:
            raise RuntimeError("Flagger executed and has not been reset. Cannot primte with new test data.")
        if not isinstance(var_params, VarParams):
            raise TypeError("Flagger can only be primed with a VarParams instance")
        if self._var_params is not None and not overwrite:
            raise RuntimeError("Flagger already primed, and overwrite is False")
        self._var_params = var_params

    def reset(
        self
    ):
        self._var_params = None
        self._executed = False

    # override if required
    def prefilter(
        self
    ):
        """
        filter reads prior to test
        """

    # modify unfrozen elements of params at runtime (or create a new params) to modify
    # TODO: AUTOMATED FILTER VERIFICATION - DEEPLY CHECK A TEST METHOD DOES NOT MODIFY ELEMENTS (if we can't restrict)
    @abstractmethod
    def test(
        self,
    ) -> FilterResult_T:   # may not modify reads
        """
        Abstract filter test method for subclass override.
        
        Each filter must define a test method via override, respecting the method signature.

        Filter test methods are for applying logic on the read data relevant to a variant, to decide whether to mark
        the variant in the FILTER field of an output VCF (or use the result in some other way).

        The test method of any subclass must return a tuple containing an object of type ReadCollection_T (Collection[AlignedSegment]
        or Mapping[Any, AlignedSegment]), and an object of type FilterResult_T (a specific subclass of FilterResult),
        where the specific types are set according to the FilterTester subclass definition. In simple terms, the test method must
        return both reads, and the result of testing. The reason for this as follows - since some filters are expected
        to mutate/drop reads from analysis as to exclude them for analysis in downstream filters, this enforced return type
        ensures that the returned value of the test method can be used in largely the same way regardless of whether the
        specific filter mutates the input reads. If a specific test does not need to drop reads from downstream analysis,
        simply return the input reads unmodified.

        ReadCollection_T is broad over both Collections and Mappings - it is expected that some filters don't need to be
        aware of which sample (in a multisample VCF) the reads came from, in which case you might make the subclass of
        FilterTester specific to list[AlignedSegment], and that some filters do need to be aware of which sample the reads
        came from, in which case you might make the subclass of FilterTester specifc to dict[str, AlignedSegment]
        """

    def __call__(
        self,
        var_params: VarParams_T | None = None,
        *args: str
    ):
        """
        run prefilter, test variant, and return result
        """
        # intentionally hidden dev options
        force = True if "force" in args else False
        overwrite = True if "overwrite" in args else False
        execute_then_reset = True if "execute_then_reset" else False

        if self._executed and not force:
            raise RuntimeError("Flagger executed and has not been reset! Cannot run flagger.")
        if var_params is not None:
            self.prime(var_params, overwrite)
        else:
            try:
                self.var_params
            except:
                raise RuntimeError("var_params has not been set! Cannot run flagger.")
        self.prefilter()
        result = self.test()
        self._executed = True
        if execute_then_reset == True:
            self.reset()  # disengage
        return result
