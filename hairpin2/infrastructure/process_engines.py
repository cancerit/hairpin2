import enum
from abc import abstractmethod
from collections.abc import Sequence
from typing import Callable, Protocol, TypeVar, runtime_checkable

from hairpin2.infrastructure.process_params import FixedParams, RunParams
from hairpin2.infrastructure.structures import FlagResult

RunParams_T = TypeVar("RunParams_T", bound=RunParams)
FixedParams_T = TypeVar("FixedParams_T", bound=FixedParams, covariant=True)
EngineResult_T = TypeVar("EngineResult_T", bound=FlagResult, covariant=True)
OptFixedParams_T = TypeVar("OptFixedParams_T", bound=FixedParams | None, covariant=True)
OptEngineResult_T = TypeVar("OptEngineResult_T", bound=FlagResult | None, covariant=True)


@runtime_checkable
class ProcessEngineProtocol[R: RunParams, O: FlagResult | None](Protocol):
    @abstractmethod
    def run_process(  # pyright: ignore[reportInvalidAbstractMethod]
        self, run_params: R
    ) -> O:
        raise NotImplementedError


# TODO: may also want to modify variants, filter variants...
# more engines to come


class ReadTaggerEngine[R: RunParams, F: FixedParams | None]:
    _adds_marks: set[str]
    _mod_func: Callable[[R, FixedParams], None] | Callable[[R, None], None]
    _read_proc_params: F

    def __init__(
        self,
        params: F,
        mod_func: Callable[[R, FixedParams_T], None] | Callable[[R, None], None],
        adds_marks: Sequence[str],
    ) -> None:
        self._read_proc_params = params
        self._mod_func = mod_func  # pyright: ignore[reportAttributeAccessIssue]
        self._adds_marks = set(adds_marks)

    def run_process(self, run_params: R) -> None:
        self._mod_func(run_params, self.tagger_params)  # pyright: ignore[reportArgumentType]

    @property
    def tagger_params(self) -> F:
        return self._read_proc_params


class VariantFlaggerEngine:
    _result_type: type[FlagResult]
    _flagger_func: Callable[[RunParams, FixedParams], FlagResult]
    _flagger_params: FixedParams

    def __init__(
        self,
        params: FixedParams_T,
        flagger_func: Callable[[RunParams_T, FixedParams_T], EngineResult_T],
        result_type: type[EngineResult_T],
    ) -> None:
        self._flagger_func = flagger_func  # pyright: ignore[reportAttributeAccessIssue]
        self._result_type = result_type

        self._flagger_params = params

    def run_process(self, run_params: RunParams) -> FlagResult:
        res = self._flagger_func(run_params, self.flagger_params)
        if not isinstance(res, self._result_type):
            raise RuntimeError("Flagger returned wrong result type - engine misconfigured")
        return res

    @property
    def flagger_params(self) -> FixedParams:
        return self._flagger_params


@enum.unique
class ProcessKindEnum(enum.Enum):
    TAGGER = ReadTaggerEngine
    FLAGGER = VariantFlaggerEngine
