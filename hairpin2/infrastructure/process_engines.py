# pyright: reportImplicitStringConcatenation=false

from abc import abstractmethod
from collections.abc import Mapping, Sequence
import enum
from typing import Callable, Protocol, TypeVar, Any, cast, override, runtime_checkable

from hairpin2.infrastructure.process_params import FixedParams, RunParams
from hairpin2.infrastructure.structures import FlagResult


# TODO: docstrings completely outdated
# TODO: consider using decorator to inspect sigs of scientific functions and wrap them
# allowing for scientific logic to be near independent of knowing anything
# about this package. Oh if possible this would be great - could extract args and typehints
# and use that to request fixed params, and even runtime params, without
# definition of a param dataclass. THIS IS NICE
# Also, since the type reqs for doing this is simply that an object is callable,
# for any stateful approaches across multiple records/positions/bundles of reads
# just use a class def with a __call__ method for advanced users
#


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


RunParams_contraT = TypeVar("RunParams_contraT", bound=RunParams, contravariant=True)
FixedParams_T = TypeVar("FixedParams_T", bound=FixedParams)
FlagResult_T = TypeVar("FlagResult_T", bound=FlagResult)
EngineResult_T = TypeVar("EngineResult_T", None, FlagResult, covariant=True)


@runtime_checkable
class ProcessEngineProtocol(Protocol[EngineResult_T]):
    @abstractmethod
    def run_process(  # pyright: ignore[reportInvalidAbstractMethod]
        self, run_params: RunParams
    ) -> EngineResult_T:
        raise NotImplementedError


# TODO: may also want to modify variants, filter variants...
# more engines to come


# A bit hacky, but functional for now
class ReadTaggerEngine(ProcessEngineProtocol[None]):
    _adds_marks: set[str]
    _param_class: type[FixedParams] | None
    _parameterised_modifier: Callable[[RunParams, FixedParams], None] | None
    _parameterless_modifier: Callable[[RunParams], None] | None
    _read_proc_params: FixedParams | None

    def __init__(
        self,
        params: Mapping[str, Any],
        param_class: type[FixedParams] | None,
        mod_func: Callable[[RunParams_contraT, FixedParams_T], None]
        | Callable[[RunParams_contraT], None],
        adds_marks: Sequence[str],
    ) -> None:
        self._param_class = param_class
        if param_class is None:
            self._parameterless_modifier = cast(Callable[[RunParams], None], mod_func)
            self._parameterised_modifier = None
        else:
            self._parameterless_modifier = None
            self._parameterised_modifier = cast(Callable[[RunParams, FixedParams], None], mod_func)
        self._adds_marks = set(adds_marks)

        if self._param_class is not None:
            requisite = {fd for fd in self._param_class.model_fields}
            kv = {k: v for k, v in params.items() if k in requisite}
            self._read_proc_params = self._param_class(**kv)
        else:
            self._read_proc_params = None

    @override
    def run_process(self, run_params: RunParams) -> None:
        if self._parameterised_modifier is not None:
            self._parameterised_modifier(run_params, cast(FixedParams, self.tagger_params))
        elif self._parameterless_modifier is not None:
            self._parameterless_modifier(run_params)

    @property
    def tagger_params(self) -> FixedParams | None:
        return self._read_proc_params


class VariantFlaggerEngine(ProcessEngineProtocol[FlagResult_T]):
    _param_class: type[FixedParams]
    _result_type: type[FlagResult]
    _flagger_func: Callable[[RunParams, FixedParams], FlagResult]
    _flagger_params: FixedParams

    def __init__(
        self,
        params: Mapping[str, Any],
        param_class: type[FixedParams],
        flagger_func: Callable[[RunParams_contraT, FixedParams_T], FlagResult_T],
        result_type: type[FlagResult],
    ) -> None:
        self._param_class = param_class
        self._flagger_func = cast(
            Callable[[RunParams, FixedParams], FlagResult_T], flagger_func
        )  # TODO: figure out if these casts are sensible
        self._result_type = result_type

        requisite = {fd for fd in self._param_class.model_fields}
        kv = {k: v for k, v in params.items() if k in requisite}
        self._flagger_params = self._param_class(**kv)

    @override
    def run_process(self, run_params: RunParams) -> FlagResult:
        res = self._flagger_func(run_params, self.flagger_params)
        if not isinstance(res, self._result_type):
            raise RuntimeError("Flagger returned wrong result type - engine misconfigured")
        return res

    @property
    def flagger_params(self) -> FixedParams:
        return self._flagger_params


# END SECTION: Behaviour Mixins -------------------------


@enum.unique
class ProcessTypeEnum(enum.Enum):
    TAGGER = ReadTaggerEngine
    FLAGGER = VariantFlaggerEngine
