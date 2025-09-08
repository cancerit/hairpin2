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

# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportIncompatibleVariableOverride=false
# pyright: reportUnnecessaryIsInstance=false
# pyright: reportImplicitStringConcatenation=false
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
from types import new_class
from pydantic import BaseModel, ConfigDict, field_validator
from typing import Callable, ClassVar, TypeVar, Any, cast, overload, override
from pysam import AlignedSegment, VariantRecord
from hairpin2.abstractions.structures import ReadView


# TODO/BUG: docstrings completely outdated
# TODO: move data structures into abstractions.structures




# SECTION: Params --------------------------------
# Users define the params required for a process
# using the nice pydantic dataclass-like syntax
# and then they can use those in their process functions
# with dot syntax. nice. Also free pydantic validation

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
class FixedParams(_Params):
    pass


# some mandatory fields
class RunParams(_Params):
    record: VariantRecord
    reads: ReadView


# END SECTION: Params --------------------------------


# SECTION: Behaviour Mixins --------------------------------
# Mixin class/es per behaviour, and a decorator (factory)
# to configure the class and inject it into user-defined subclasses
# of ReadAwareProcess. The point is to let the user
# declare read aware processes in a modular way, without
# having to deeply understand generics, inheritance, overrides, etc.
# and allowing methods to be written as free functions
# for easy export and testing. This approach also allows
# for pretty effective leveraging of the type checker compared to
# e.g. __init_subclass__.
# TODO: decorators could have better nomenclature, so they look less like regular funcs


# TODO: use elsewhere, complete typing
RunParams_contraT = TypeVar("RunParams_contraT", bound=RunParams, contravariant=True)
FixedParams_T = TypeVar("FixedParams_T", bound=FixedParams)
C = TypeVar("C", bound=type)


# TODO: extend, subsume into other classes
class _RequireReadProperties:
    RequireTags: ClassVar[tuple[str, ...]]
    ExcludeTags: ClassVar[tuple[str, ...]]
    RequireFields: ClassVar[tuple[str, ...]]

    def require_properties_check(
        self,
        run_params: RunParams
    ):
        """
        filter reads prior to test
        """
        # TODO: global on-fail options: record/split offending reads to file, record/split variant to file, ignore, warn, fail
        # TODO: check primed first!

        # TODO: doc - you don't need to use RequireReadProperties, but it injects nice behaviour for you
        filtered: dict[Any, list[AlignedSegment]] = {}
        for sample, reads in run_params.reads.items():
            passed_reads: list[AlignedSegment] = []
            for read in reads:
                rpass = True
                if self.RequireTags and not all(read.has_tag(tag) for tag in self.RequireTags):
                    rpass = False
                if self.ExcludeTags and any(read.has_tag(tag) for tag in self.ExcludeTags):
                    rpass = False
                if self.RequireFields and any(getattr(read, field, None) is None for field in self.RequireFields):
                    rpass = False
                if rpass:
                    passed_reads.append(read)

            filtered[sample] = passed_reads
        return ReadView(filtered)


# decorator
def require_read_properties(
    *,
    require_tags: Iterable[str] | None = None,
    exclude_tags: Iterable[str] | None = None,
    require_fields: Iterable[str] | None = None
):
    def deco[C: type](cls: C) -> C:
        if issubclass(cls, _RequireReadProperties):
            raise TypeError(
                f"{cls.__name__} is already of type {_RequireReadProperties}, "
                "perhaps you used the require_read_properties decorator twice?"
            )
        bases = (cls, _RequireReadProperties)
        raf_subclass_kwds  = getattr(cls, '__decl_kwargs__', {}) 

        def body(ns: dict[str, Any]):
            ns["__module__"] = cls.__module__
            ns["__doc__"] = getattr(cls, "__doc__", None)
            ns["__qualname__"] = getattr(cls, "__qualname__", cls.__name__)

        new_cls = cast(C, new_class(cls.__name__, bases, raf_subclass_kwds, body))
        setattr(new_cls, "RequireTags", tuple(require_tags) if require_tags else tuple())
        setattr(new_cls, "ExcludeTags", tuple(exclude_tags) if exclude_tags else tuple())
        setattr(new_cls, "RequireFields", tuple(require_fields) if require_fields else tuple())
        return new_cls
    return deco


# TODO: may also want to modify variants, filter variants...
# more mixins to come


class _PrefilterProcess:
    PrefilterParamClass: ClassVar[type[FixedParams]]
    _ReadEvaluator: ClassVar[Callable[[AlignedSegment, VariantRecord, FixedParams], bool]]
    _prefilter_params: FixedParams | None = None

    def load_prefilter(self, *, params: Mapping[str, Any]) -> None:
        requisite = {fd for fd in self.PrefilterParamClass.model_fields}
        kv = {k: v for k, v in params.items() if k in requisite}
        self._prefilter_params = self.PrefilterParamClass(**kv)

    # “sensible default” implementation; override if needed - # BUG: well you can't if it's _private! (ok you can but you shouldn't)
    def filter_reads(self, run_params: RunParams) -> Mapping[Any, Sequence[AlignedSegment]]:
        if self.prefilter_params is None:
            raise RuntimeError("Prefilter not loaded; call load_prefilter() first")
        out: dict[Any, list[AlignedSegment]] = {}
        for sample_key, reads in run_params.reads.items():
            passed: list[AlignedSegment] = []
            for read in reads:
                if type(self)._ReadEvaluator(read, run_params.record, self.prefilter_params):
                    passed.append(read)
            out[sample_key] = passed
        return out

    @property
    def prefilter_params(self) -> FixedParams | None:
        return self._prefilter_params


# decorator
def prefilter[FixedParams_T: FixedParams](
    *,
    prefilter_param_class: type[FixedParams_T],
    read_evaluator_func: Callable[[AlignedSegment, VariantRecord, FixedParams_T], bool],
):
    def deco[C: type](cls: C) -> C:
        if issubclass(cls, _PrefilterProcess):
            raise TypeError(
                f"{cls.__name__} is already of type {_PrefilterProcess}, "
                "perhaps you used the prefilter decorator twice?"
            )

        raf_subclass_kwds  = getattr(cls, '__decl_kwargs__', {}) 
        bases = (cls, _PrefilterProcess)

        def body(ns: dict[str, Any]):
            ns["__module__"] = cls.__module__
            ns["__doc__"] = getattr(cls, "__doc__", None)
            ns["__qualname__"] = getattr(cls, "__qualname__", cls.__name__)

        new_cls = cast(C, new_class(cls.__name__, bases, raf_subclass_kwds, body))

        setattr(new_cls, "PrefilterParamClass", prefilter_param_class)
        setattr(new_cls, "_ReadEvaluator", read_evaluator_func)
        return new_cls
    return deco


# TODO: all mixins should probably allow for both fixed params and no fixed params subtypes
class _AbstractReadTaggerProcess(ABC):
    AddsTag: ClassVar[str]
    ReadTaggerParamClass: ClassVar[type[FixedParams] | None]
    _ReadModifier: ClassVar[Callable[[RunParams, FixedParams], None] | Callable[[RunParams], None]]
    _read_proc_params: FixedParams | None = None
    _loaded: int = False

    def load_tagger(self, *, params: Mapping[str, Any]) -> None:
        if self.ReadTaggerParamClass is not None:
            requisite = {fd for fd in self.ReadTaggerParamClass.model_fields}
            kv = {k: v for k, v in params.items() if k in requisite}
            self._read_proc_params = self.ReadTaggerParamClass(**kv)
        else:
            self._read_proc_params = None
        self._loaded = True

    # TODO: it would be safer for the tagger to hand back a bool and this mixin handle the tagging
    @abstractmethod
    def modify_reads(self, run_params: RunParams) -> None:
        ...

    @property
    @abstractmethod
    def tagger_params(self) -> FixedParams | None:
        ...


class _FixedParamsReadTaggerProcess(_AbstractReadTaggerProcess):
    AddsTag: ClassVar[str]
    ReadTaggerParamClass: ClassVar[type[FixedParams]]
    _ReadModifier: ClassVar[Callable[[RunParams, FixedParams], None]]
    _read_proc_params: FixedParams
    _loaded: int = False

    @override
    def modify_reads(self, run_params: RunParams) -> None:
        if self.tagger_params is None:
            raise RuntimeError("Read tagger not loaded; call load_tagger() first")
        type(self)._ReadModifier(run_params, self.tagger_params)

    @property
    @override
    def tagger_params(self) -> FixedParams | None:
        return self._read_proc_params


class _NoParamsReadTaggerProcess(_AbstractReadTaggerProcess):
    AddsTag: ClassVar[str]
    ReadTaggerParamClass: ClassVar[None] = None
    _ReadModifier: ClassVar[Callable[[RunParams], None]]
    _read_proc_params: None = None
    _loaded: int = False

    @override
    def modify_reads(self, run_params: RunParams) -> None:
        if not self._loaded:
            raise RuntimeError("Read tagger not loaded; call load_tagger() first")
        type(self)._ReadModifier(run_params)

    @property
    @override
    def tagger_params(self) -> None:
        return self._read_proc_params


# decorator
@overload
def read_tagger(
    *,
    tagger_param_class: None,
    read_modifier_func: Callable[[RunParams_contraT], None],
    adds_tag: str
) -> Callable[[C], C]: ...
@overload
def read_tagger(
    *,
    tagger_param_class: type[FixedParams_T],
    read_modifier_func: Callable[[RunParams_contraT, FixedParams_T], None],
    adds_tag: str
) -> Callable[[C], C]: ...

def read_tagger(
    *,
    tagger_param_class: type[FixedParams_T] | None,
    read_modifier_func: Callable[[RunParams_contraT, FixedParams_T], None] | Callable[[RunParams_contraT], None],
    adds_tag: str
) -> Callable[[C], C]:
    def deco(cls: C) -> C:
        if issubclass(cls, _AbstractReadTaggerProcess):
            raise TypeError(
                f"{cls.__name__} is already of type {_AbstractReadTaggerProcess}, "
                "perhaps you used the read_tagger decorator twice, or inherited as well as decorated?"
            )

        if tagger_param_class is None:
            bases = (cls, _NoParamsReadTaggerProcess)
        else:
            bases = (cls, _FixedParamsReadTaggerProcess)
        raf_subclass_kwds  = getattr(cls, '__decl_kwargs__', {}) 

        def body(ns: dict[str, Any]):
            ns["__module__"] = cls.__module__  # so user class comes from the same module
            ns["__doc__"] = getattr(cls, "__doc__", None)  # ostensibly good pratice (shrug)
            ns["__qualname__"] = getattr(cls, "__qualname__", cls.__name__)

        new_cls = cast(C, new_class(cls.__name__, bases, raf_subclass_kwds, body))

        setattr(new_cls, "ReadTaggerParamClass", tagger_param_class)
        setattr(new_cls, "_ReadModifier", read_modifier_func)
        setattr(new_cls, "AddsTag", adds_tag)
        return new_cls
    return deco


# NOTE/TODO: FlagResult is the part of the abstractflaggers model about which I am most skeptical
# since it uses init_subclass over a decorator, making it different
# and because it requires the user to override an abstract method
class FlagResult(BaseModel, ABC):  # pyright: ignore[reportUnsafeMultipleInheritance]
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
    """
    FlagName: ClassVar[str]
    AllowedCodes: ClassVar[tuple[int, ...] | tuple[str, ...] | None]
    CodeType: ClassVar[type]
    flag: bool | None
    code: int | str | None = None

    model_config: ConfigDict = ConfigDict(
        strict=True,
        frozen=True,
        arbitrary_types_allowed=True
    )

    @field_validator('code', mode='after')
    @classmethod
    def validate_code(
        cls,
        value: int | str | None
    ):
        if cls.AllowedCodes is None and value is not None:
            raise AttributeError(f"The class definition for FlagResult {cls.FlagName!r} does not specify result codes, so you cannot set a result code")
        elif cls.AllowedCodes is not None and value is None:
            raise AttributeError(f"The class definition for FlagResult {cls.FlagName!r} does specifies result codes, so you must set a result code")
        elif not isinstance(value, cls.CodeType):
            raise TypeError(f"Attempting to set result code for FlagResult {cls.FlagName!r} with value {value} of type {type(value)}, but this FlagResult specifies codes must be of type {cls.CodeType}")
        elif value is not None and cls.AllowedCodes is not None and value not in cls.AllowedCodes:
            raise ValueError(f"Attempting to set code with value {value}, but this FlagResult {cls.FlagName!r} only allows {cls.AllowedCodes}")
        return value

    def __init_subclass__(
        cls,
        *,
        flag_name: str,
        result_codes: Iterable[int] | Iterable[str] | None,
        **kwargs: Any
    ) -> None:
        # if 'abstractflaggers' in cast(str, getattr(cls, "__module__", None)):
        #     # This is an intermediate pydantic class,
        #     # not a user result class
        #     # skip or everything breaks
        #     return

        cls.FlagName = flag_name

        match result_codes:
            case None:
                cls.AllowedCodes = None
            case [*items] if all(isinstance(el, int) for el in items):
                cls.AllowedCodes = cast(tuple[int, ...], tuple(result_codes))
                cls.CodeType = int
            case [*items] if all(isinstance(el, str) for el in items):
                cls.AllowedCodes = cast(tuple[str, ...], tuple(result_codes))
                cls.CodeType = str
            case _:
                raise TypeError(f"allowed_codes for FlagResult {cls.FlagName} provided must be an iterable of all int, or all str; or None")
        
        super().__init_subclass__(**kwargs)

    @abstractmethod
    def getinfo(self) -> str:
        """
        Return basic filter info in a string formatted for use in the VCF INFO field - "<flag>|<code>".

        Each filter must return INFO as it should be formatted for the VCF INFO field, or None if not applicable.
        Subclasses must override this method to return more specific info.
        """


FlagResult_T = TypeVar('FlagResult_T', bound=FlagResult)


class _VariantFlagProcess:
    FlagName: ClassVar[str]
    FlaggerParamClass: ClassVar[type[FixedParams]]
    FlaggerResultType: ClassVar[type[FlagResult]]
    _Flagger: ClassVar[Callable[[RunParams, FixedParams], FlagResult]]
    _flagger_params: FixedParams | None = None

    def load_flagger(self, *, params: Mapping[str, Any]) -> None:
        requisite = {fd for fd in self.FlaggerParamClass.model_fields}
        kv = {k: v for k, v in params.items() if k in requisite}
        self._flagger_params = self.FlaggerParamClass(**kv)  # type: ignore[call-arg]

    def flag_variant(self, run_params: RunParams) -> FlagResult:
        if self.flagger_params is None:
            raise RuntimeError("Flagger not loaded; call load_flagger() first")
        res = type(self)._Flagger(run_params, self.flagger_params)
        if not isinstance(res, self.FlaggerResultType):
            raise RuntimeError("Flagger returned wrong result type")
        return res

    @property
    def flagger_params(self) -> FixedParams | None:
        return self._flagger_params

# decorator
def variant_flagger(
    *,
    flag_name: str,
    flagger_param_class: type[FixedParams_T],
    flagger_func: Callable[[RunParams_contraT, FixedParams_T], FlagResult_T],
    result_type: type[FlagResult_T],
):
    def deco[C: type](cls: C) -> C:
        # avoid double-wrapping if already present
        if issubclass(cls, _VariantFlagProcess):
            raise TypeError(
                f"{cls.__name__} is already of type {_VariantFlagProcess}, "
                "perhaps you used the decorator twice, or inherited as well as decorated?"
            )
        bases = (cls, _VariantFlagProcess)
        raf_subclass_kwds  = getattr(cls, '__decl_kwargs__', {}) 

        def body(ns: dict[str, Any]):
            ns["__module__"] = cls.__module__
            ns["__doc__"] = getattr(cls, "__doc__", None)
            ns["__qualname__"] = getattr(cls, "__qualname__", cls.__name__)

        new_cls = cast(C, new_class(cls.__name__, bases, raf_subclass_kwds, body))

        # runtime sanity checks
        if getattr(result_type, "FlagName", None) != flag_name:
            raise TypeError("flag_name must match result_type.FlagName")

        setattr(new_cls, "FlagName", flag_name)
        setattr(new_cls, "FlaggerParamClass", flagger_param_class)
        setattr(new_cls, "FlaggerResultType", result_type)
        setattr(new_cls, "_Flagger", flagger_func)
        return new_cls
    return deco


# END SECTION: Behaviour Mixins -------------------------


# The key class which executes the behaviour injected via the mixins
# TODO: mixin decorators should confirm class they are decorating
# is a ReadAwareProcess
# TODO: fix call typing
# TODO: should allow no parameter init if appropriate
class ReadAwareProcess(ABC):
    # TODO: update docstring
    _ProcessName: ClassVar[str]
    __slots__ = ("_param_map", "_var_params", "_executed")  # pyright: ignore[reportUnannotatedClassAttribute]

    def __init_subclass__(
        cls,
        *,
        process_name: str,
       **kwargs: Any
    ) -> None:
        # I don't know how robust this is but it does work
        cls._ProcessName = process_name  # swallow process name
        super().__init_subclass__(**kwargs)  # pass on the rest
        kwargs.update({'process_name': process_name})
        cls.__decl_kwargs__: dict[str, Any] = dict(kwargs)   # stash for decorators

    def __init__(
        self,
        params: Mapping[str, Any],
        **kwargs: Any
    ):

        # will fail if key not present
        # TODO: consider how to handle that
        my_params = params[self._ProcessName]
        
        self._param_map: Mapping[str, Any] = my_params
        self._var_params: RunParams | None = None
        self._executed: bool = False

        # TODO/NOTE: these all do the same thing
        if isinstance(self, _PrefilterProcess):
            self.load_prefilter(params=my_params)
        if isinstance(self, _AbstractReadTaggerProcess):
            self.load_tagger(params=my_params)
        if isinstance(self, _VariantFlagProcess):
            self.load_flagger(params=my_params)

        super().__init__(**kwargs)

    @property
    def all_fixed_params(self) -> Mapping[str, Any]:
        return self._param_map

    @property
    def run_params(self) -> RunParams:
        if self._var_params is None:
            raise AttributeError("run_params not set, cannot access")
        return self._var_params

    def prime(
        self,
        var_params: RunParams,
        overwrite: bool = False
    ):
        if self._executed:
            raise RuntimeError("Process executed and has not been reset. Cannot prime with new test data.")
        if not isinstance(var_params, RunParams):
            # The following message is not stricly true, but should be
            # TODO: init subclass should tie a type[RunParams] to the class
            raise TypeError("Process can only be primed with an instance of the RunParams subclass tied to this flagger.")  # pyright: ignore[reportUnreachable]
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
        if self._var_params is None:
            raise RuntimeError("Attempt to set filtered reads from runtime params, but runtime params is not set, so what were you filtering?!?")
        if isinstance(new_reads, dict):
            new_reads = ReadView(new_reads)
        run_arg_d = self._var_params.__dict__.copy()  # must copy, or you'll alter the parent readview!!
        run_arg_d['reads'] = new_reads  # overwrite only the readview aspect
        self.prime(self._var_params.__class__(**run_arg_d), overwrite=True)

    def reset(
        self
    ):
        self._var_params = None
        self._executed = False


    # @overload
    # def __call__(
    #     self: VariantFlagProcess  # BUG: these don't work because VariantFlagProcess isn't a subclass of this base
    # ) -> FlagResult: ...
    # @overload
    # def __call__(
    #     self
    # ) -> None: ...
    
    def __call__(
        self,
        call_run_params: RunParams | None = None,
        *,
        _internal_switches: list[str] | None = None,  # hidden dev options
    ) -> FlagResult | None:
        """
        run prefilter, process reads
        """

        switches = _internal_switches or []
        force = True if "force" in switches else False
        overwrite = True if "overwrite" in switches else False
        execute_then_reset = True if "execute_then_reset" in switches else False

        flag_result: FlagResult | None = None

        if self._executed and not force:
            # TODO: use subclass name via fstring
            raise RuntimeError("Process executed, and has not been reset and loaded with new data! Cannot run process.")
        if call_run_params is not None:
            self.prime(call_run_params, overwrite)
        else:
            try:
                self.run_params
            except:
                raise RuntimeError("Process has not been primed with data! Cannot run process.")

        # might need to use name mangling to ensure no clashes
        if isinstance(self, _RequireReadProperties):
            self._set_filtered_reads(self.require_properties_check(self.run_params))
        if isinstance(self, _PrefilterProcess):
            self._set_filtered_reads(self.filter_reads(self.run_params))
        if isinstance(self, _AbstractReadTaggerProcess):
            self.modify_reads(self.run_params)
        if isinstance(self, _VariantFlagProcess):
            flag_result = self.flag_variant(self.run_params)

        self._executed = True
        if execute_then_reset == True:
            self.reset()  # disengage

        return flag_result
