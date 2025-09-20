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

# pyright: reportImplicitStringConcatenation=false

from abc import ABC
from collections.abc import Mapping
from typing import Any, Callable, ClassVar, Protocol
from collections.abc import Sequence

from hairpin2.abstractions.process_engines import EngineResult_T, ProcessEngineProtocol, ProcessTypeEnum
from hairpin2.abstractions.process_params import RunParams
from hairpin2.abstractions.structures import ExtendedRead, FlagResult, ReadView, read_has_mark


class ReadAwareProcessProtocol(Protocol[EngineResult_T]):
    ProcessNamespace: ClassVar[str]
    EngineFactory: ClassVar[Callable[[Mapping[str, Any]], ProcessEngineProtocol[Any]]]
    ProcessType: ClassVar[ProcessTypeEnum]
    AddsMarks: ClassVar[set[str] | None]

    def __init__(
        self,
        process_namespace_fixed_params: Mapping[str, Any],
        require_marks: Sequence[str],
        exclude_marks: Sequence[str],
    ): ...

    @property
    def require_marks(
        self
    ) -> set[str]: ...

    @property
    def exclude_marks(
        self
    ) -> set[str]:
        ...

    @property
    def fixed_params(self) -> Mapping[str, Any]: ...

    @property
    def run_params(self) -> RunParams: ...

    def reset(
        self
    ) -> None: ...

    def __call__(
        self,
        call_run_params: RunParams | None = None,
        *,
        _internal_switches: Sequence[str] | None = None,  # hidden dev options
    ) -> EngineResult_T: ...


class ReadAwareProcess(ABC):
    ProcessNamespace: ClassVar[str | None] = None
    EngineFactory: ClassVar[Callable[[Mapping[str, Any]], ProcessEngineProtocol[Any]] | None] = None
    ProcessType: ClassVar[ProcessTypeEnum | None] = None
    AddsMarks: ClassVar[set[str] | None] = None  # ok for now, but really I'd like this to be stored on the engine
    # TODO: update docstring
    __slots__: tuple[str, ...] = ("_param_map", "_var_params", "_executed", "_require_marks", "_exclude_marks", "_engine")

    def __init__(
        self,
        process_namespace_fixed_params: Mapping[str, Any],
        require_marks: Sequence[str],
        exclude_marks: Sequence[str],
    ):
        cls = type(self)
        if any(x is None for x in (cls.ProcessNamespace, cls.EngineFactory, cls.ProcessType)):
            raise TypeError('Process not configured! Missing EngineFactory, ProcessNamespace, or ProcessType')
        if cls.ProcessType == ProcessTypeEnum.TAGGER and cls.AddsMarks is None:
            raise TypeError('taggers must declare added marks')

        self._param_map: Mapping[str, Any] = process_namespace_fixed_params
        self._var_params: RunParams | None = None
        self._executed: bool = False

        req_set = set(require_marks)
        exc_set = set(exclude_marks)
        if req_set & exc_set:
            raise ValueError("require_marks and exclude_marks overlap!")
        self._require_marks: set[str] = req_set
        self._exclude_marks: set[str] = exc_set

        # init engine
        assert cls.EngineFactory is not None  # type checker..., and refactor safeguard
        self._engine: ProcessEngineProtocol[Any] = cls.EngineFactory(self.fixed_params)

        if not isinstance(self._engine, ProcessEngineProtocol):
            raise RuntimeError('Engine does not appear to satisfy necessay protocol, Process misconfigured')

    def __init_subclass__(
        cls,
        process_namespace: str | None = None,
        engine_factory: Callable[[Mapping[str, Any]], ProcessEngineProtocol[Any]] | None = None,
        process_type: ProcessTypeEnum | None = None,
        adds_marks: set[str] | None = None
    ):
        cls.ProcessNamespace = process_namespace
        cls.EngineFactory = engine_factory
        cls.ProcessType = process_type
        cls.AddsMarks = adds_marks

    # since marks are based only on presence
    # the lack of an exclude mark is enough to pass a check
    # this makes it easy to enable and disable tagger processes
    # while leaving it up to the user which way round they want to use a tag
    # e.g. if you want to be able to easily turn something on/off
    # make it an exclude mark like LOW-QUAL
    # conversely if you want to be strict make a require mark like HIGH-QUAL
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
        filtered: dict[Any, list[ExtendedRead]] = {}
        for sample, reads in run_params.reads.items():
            passed_reads: list[ExtendedRead] = []
            for read in reads:
                rpass = True
                if self.require_marks and not all(read_has_mark(read, mark) for mark in self.require_marks):
                    rpass = False
                if self.exclude_marks and any(read_has_mark(read, mark) for mark in self.exclude_marks):
                    rpass = False
                if rpass:
                    passed_reads.append(read)

            filtered[sample] = passed_reads
        return ReadView(filtered)

    @property
    def require_marks(
        self
    ):
        return self._require_marks

    @property
    def exclude_marks(
        self
    ):
        return self._exclude_marks

    @property
    def fixed_params(self) -> Mapping[str, Any]:
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
        new_reads: Mapping[Any, Sequence[ExtendedRead]] | ReadView[ExtendedRead] | None
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
    
    def __call__(
        self,
        call_run_params: RunParams | None = None,
        *,
        _internal_switches: Sequence[str] | None = None,  # hidden dev options
    ) -> FlagResult | None:
        """
        run prefilter, process reads
        """

        switches = _internal_switches or []
        force = True if "force" in switches else False
        overwrite = True if "overwrite" in switches else False
        execute_then_reset = True if "execute_then_reset" in switches else False

        if self._executed and not force:
            raise RuntimeError(f"Process {type(self).__name__} with namespace {type(self).ProcessNamespace} has been executed, and has not been reset and loaded with new data! Cannot run process.")
        if call_run_params is not None:
            self.prime(call_run_params, overwrite)
        else:
            try:
                self.run_params
            except:
                raise RuntimeError(f"Process {type(self).__name__} with namespace {type(self).ProcessNamespace} has not been primed with data! Cannot run process.")

        if self.require_marks or self.exclude_marks:
            self._set_filtered_reads(self.require_properties_check(self.run_params))

        # TODO: add back prefilter functionality

        ret = self._engine.run_process(self.run_params)

        self._executed = True
        if execute_then_reset == True:
            self.reset()  # disengage

        return ret
