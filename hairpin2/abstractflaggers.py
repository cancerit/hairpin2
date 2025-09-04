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
# pyright: reportUnnecessaryIsInstance=false
"""
A set of Parent Abstract Base Classes that together completely describe the expected implementation of a scientific test
for recording FILTER entry for a variant record in a VCF

Heavily leverages static type checking in addition to runtime guarantees. As such, implementation of a new filter is quite verbose.
The payoff for that verbosity is:
    - complete static and runtime checking of filter validity, which is helpful when describing many filters, and complex tests
    - strong guarantees about the implementation of any given filter without needing to know about the underlying logic, making filters very easy to use once defined
"""
from abc import ABC, abstractmethod
from collections.abc import Iterable, Mapping, Sequence
from pydantic import BaseModel, ConfigDict
from typing import Generic, TypeVar, Any, cast, override, ClassVar
from pysam import AlignedSegment, VariantRecord
from enum import IntEnum
from hairpin2.structures import ReadView


# TODO: test instantiating subclass from json


### READ FilterTester First and work backwards! ###

CodeEnum_T = TypeVar("CodeEnum_T", bound=IntEnum, covariant=True)  # TODO: ideally minimise extra imports for definer
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
    Name: ClassVar[str]  # for VCF FILTER field  # TODO: name via __init_subclass__ methods instead
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


class _Params(BaseModel):
    model_config: ConfigDict = ConfigDict(
        strict=True,
        frozen=True,
        arbitrary_types_allowed=True
    )

    '''
    parent dataclass to be be inherited from to store specific fixed parameters for a particular subclass of FilterTester,
    or in other words for a particular filtering test. Using subclasses of this class for the fixed parameters provides
    type-safety and a consistent interface for implementing filters
    '''
    pass


# inherit and set defaults in subclasses
# exclude_flag: int
# require_flag: int
# min_avg_clip_quality: int   # in subclass
# so basically an implementation (like hp2) entirely defines prefiltering beyond checking field and tag presence
# and those subclasses define the json config schema
# myobj.model_validate(json_str)
class PrefilterParams(_Params):
    pass


class FixedParams(_Params):
    pass


class RunParams(_Params):
    record: VariantRecord
    reads: ReadView


PrefilterParams_T = TypeVar("PrefilterParams_T", bound=PrefilterParams)
FixedParams_T = TypeVar("FixedParams_T", bound=FixedParams)
RunParams_T = TypeVar("RunParams_T", bound=RunParams)
FilterResult_T = TypeVar("FilterResult_T", bound=FilterResult[IntEnum], covariant=True)  # covariant such that a test method that returns a subtype of FilterResult[IntEnum] is accepted where FilterResult[IntEnum] (or FilterResult_T) is expected


class RequireReadProperties:
    def __init_subclass__(
        cls,
        *,
        require_tags: Iterable[str],
        exclude_tags: Iterable[str],
        require_fields: Iterable[str],
        **kwargs: Any
    ):
        super().__init_subclass__(**kwargs)

        if cls is RequireReadProperties:
            return

        cls.require_tags: tuple[str, ...] = tuple(require_tags)
        cls.exclude_tags: tuple[str, ...] = tuple(exclude_tags)
        cls.require_fields: tuple[str, ...] = tuple(require_fields)


# TODO: inject it all. subclass unnecssary, is processor or tester depending on mixins only, and prefilter
# way less for end user to define, easier to share behaviour and not duplicate code
# way easier to test as can just test the callable without the machinery
class _BaseExecutor(Generic[PrefilterParams_T, FixedParams_T, RunParams_T], ABC):
    # TODO: update docstring
    __slots__ = ("_prefilter_params", "_fixed_params", "_var_params", "_executed")

    def __init__(  # pyright: ignore[reportMissingSuperCall]
        self,
        prefilter_params: PrefilterParams_T,
        fixed_params: FixedParams_T
    ):
        self._prefilter_params: PrefilterParams_T = prefilter_params
        self._fixed_params: FixedParams_T = fixed_params
        self._var_params: RunParams_T | None = None
        self._executed: bool = False

    @property
    def prefilter_params(self) -> PrefilterParams_T:
        return self._prefilter_params

    @property
    def fixed_params(self) -> FixedParams_T:
        return self._fixed_params

    @property
    def run_params(self) -> RunParams_T:
        if self._var_params is None:
            raise AttributeError("var_params not set, cannot access")
        return self._var_params

    def prime(
        self,
        var_params: RunParams_T,
        overwrite: bool = False
    ):
        if self._executed:
            raise RuntimeError("Flagger executed and has not been reset. Cannot primte with new test data.")
        if not isinstance(var_params, RunParams):
            raise TypeError("Flagger can only be primed with an instance of the RunParams subclass tied to this flagger.")  # pyright: ignore[reportUnreachable]
        if self._var_params is not None and not overwrite:
            raise RuntimeError("Flagger already primed, and overwrite is False")
        self._var_params = var_params

    # to support frozen run params
    # provide a sneaky way to rebuild post read filtering
    # without burdening the user
    def _set_filtered_reads(
        self,
        new_reads: Mapping[Any, Sequence[AlignedSegment]] | ReadView | None
    ):
        if isinstance(new_reads, dict):
            new_reads = ReadView(new_reads)
        run_arg_d = self._var_params.__dict__
        run_arg_d['reads'] = new_reads
        self.prime(cast(RunParams_T, self._var_params.__class__(**run_arg_d)), overwrite=True)

    def reset(
        self
    ):
        self._var_params = None
        self._executed = False

    def _internal_prefilter(
        self
    ):
        """
        filter reads prior to test
        """
        # TODO: global on-fail options: record/split offending reads to file, record/split variant to file, ignore, warn, fail
        # TODO: check primed first!

        # TODO: doc - you don't need to use RequireReadProperties, but it injects nice behaviour for you
        # NOTE: might be nice to inject prefiltering this way? consider. but not necessary for this release
        if isinstance(self, RequireReadProperties):
            filtered: dict[Any, list[AlignedSegment]] = {}
            for sample, reads in self.run_params.reads.items():
                passed_reads: list[AlignedSegment] = []
                for read in reads:
                    rpass = True
                    if not all(read.has_tag(tag) for tag in self.require_tags):
                        rpass = False
                    if any(read.has_tag(tag) for tag in self.exclude_tags):
                        rpass = False
                    if any(getattr(read, field, None) is None for field in self.require_fields):
                        rpass = False
                    if rpass:
                        passed_reads.append(read)

                filtered[sample] = passed_reads

            self._set_filtered_reads(ReadView(filtered))

    # user override if required
    def prefilter(
        self
    ) -> Mapping[Any, Sequence[AlignedSegment]]:
        # must check var params set
        if self._var_params is None:
            raise RuntimeError("Attempt made to prefilter without data, failed")
        return self._var_params.reads


class ReadPreprocessor(
    _BaseExecutor[PrefilterParams_T, FixedParams_T, RunParams_T],
    Generic[PrefilterParams_T, FixedParams_T, RunParams_T],
    ABC
):

    @abstractmethod
    def process(
        self
    ) -> None:  # must modify reads in place, hence None return
        """
        additive read processing - modify or add to pysam AlignedSegment objects.
        """

    def __call__(
        self,
        var_params: RunParams_T | None = None,
        *,
        _internal_switches: list[str] | None = None,  # hidden dev options
    ):
        """
        run prefilter, process reads
        """

        switches = _internal_switches or []
        force = True if "force" in switches else False
        overwrite = True if "overwrite" in switches else False
        execute_then_reset = True if "execute_then_reset" in switches else False

        if self._executed and not force:
            # TODO: use subclass name via fstring
            raise RuntimeError("Preprocessor executed, and has not been reset and loaded with new data! Cannot run preprocessor.")
        if var_params is not None:
            self.prime(var_params, overwrite)
        else:
            try:
                self.run_params
            except:
                raise RuntimeError("var_params has not been set! Cannot run preprocessor.")
        self._internal_prefilter()
        user_filtered_reads = self.prefilter()
        self._set_filtered_reads(user_filtered_reads)  # this is wasteful if prefiltering isn't overridden
        self.process()
        self._executed: bool = True
        if execute_then_reset == True:
            self.reset()  # disengage


class Flagger(
    _BaseExecutor[PrefilterParams_T, FixedParams_T, RunParams_T],
    Generic[PrefilterParams_T, FixedParams_T, RunParams_T, FilterResult_T],
    ABC
):
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

    # modify unfrozen elements of params at runtime (or create a new params) to modify
    # TODO: AUTOMATED FILTER VERIFICATION - DEEPLY CHECK A TEST METHOD DOES NOT MODIFY ELEMENTS (if we can't restrict)
    @abstractmethod
    def test(
        self,
    ) -> FilterResult_T:  # should not modify underlying reads, but not currently enforced
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
        """

    def __call__(
        self,
        var_params: RunParams_T | None = None,
        *,
        _internal_switches: list[str] | None = None,  # hidden dev options
    ):
        """
        run prefilter, test variant, and return result
        """

        switches = _internal_switches or []
        force = True if "force" in switches else False
        overwrite = True if "overwrite" in switches else False
        execute_then_reset = True if "execute_then_reset" in switches else False

        if self._executed and not force:
            raise RuntimeError("Flagger executed, and has not been reset and loaded with new data! Cannot run flagger.")
        if var_params is not None:
            self.prime(var_params, overwrite)
        else:
            try:
                self.run_params
            except:
                raise RuntimeError("var_params has not been set! Cannot run flagger.")
        self._internal_prefilter()
        user_filtered_reads = self.prefilter()
        self._set_filtered_reads(user_filtered_reads)  # this is wasteful if prefiltering isn't overridden
        result = self.test()
        self._executed: bool = True
        if execute_then_reset == True:
            self.reset()  # disengage
        return result
