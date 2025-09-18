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

from abc import abstractmethod
from collections.abc import Mapping, Sequence
from enum import Enum
from typing import Callable, Protocol, TypeVar, Any, cast, override, runtime_checkable

from hairpin2.abstractions.process_params import FixedParams, RunParams
from hairpin2.abstractions.structures import FlagResult
# from inspect import BoundArguments, Parameter, Signature, signature as inspect_sig


# TODO/BUG: docstrings completely outdated
# TODO: move data structures into abstractions.structures
# TODO: consider using wrapt more
# TODO: consider using decorator to inspect sigs of scientific functions and wrap them
# allowing for scientific logic to be near independent of knowing anything
# about this package. Oh if possible this would be great - could extract args and typehints
# and use that to request fixed params, and even runtime params, without
# definition of a param dataclass. THIS IS NICE
# Also, since the type reqs for doing this is simply that an object is callable,
# for any stateful approaches across multiple records/positions/bundles of reads
# just use a class def with a __call__ method for advanced users 
# 
# NOTE: mixin approach allows overrides in the user class definition of mixin methods
# without my injected methods stomping on them
# can still duck type regardless if desired
# TODO: use init_subclass in mixins, and pass decorator args into those as keywords
# this will allow definition time validation and guard against misspellings and stuff
# and don't set any other vars outside of exec body func
# TODO: decorators from deco factories should return protocols for the thing they're returning - won't obscure user typing, and will add type hints
# for injected behaviour - lies, doesn't work at definition time (but does at runtime, which I don't know if I want or care about)
# TODO: focus on structures lol don't get too carried away here make sure I'm doing stuff for the aim - make this easy for scientists!


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
FlagResult_T = TypeVar('FlagResult_T', bound=FlagResult)


@runtime_checkable
class ProcessEngineProtocol(Protocol):

    @abstractmethod
    def run_process(  # pyright: ignore[reportInvalidAbstractMethod]
        self,
        run_params: RunParams
    ) -> Any:
        raise NotImplementedError




# TODO: may also want to modify variants, filter variants...
# more mixins to come


# A bit hacky, but functional for now
class ReadTaggerEngine(
    ProcessEngineProtocol
):
    _adds_marks: set[str]
    _param_class: type[FixedParams] | None
    _parameterised_modifier: Callable[[RunParams, FixedParams], None] | None
    _parameterless_modifier: Callable[[RunParams], None] | None
    _read_proc_params: FixedParams | None

    def __init__(
        self,
        params: Mapping[str, Any],
        param_class: type[FixedParams] | None,
        mod_func: Callable[[RunParams_contraT, FixedParams_T], None] | Callable[[RunParams_contraT], None],
        adds_marks: Sequence[str],
    ) -> None:
        self._param_class = param_class
        if param_class is None:
            self._parameterless_modifier = cast(Callable[[RunParams], None], mod_func)
        else:
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


class VariantFlaggerEngine(
    ProcessEngineProtocol
):
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
        self._flagger_func = cast(Callable[[RunParams, FixedParams], FlagResult_T], flagger_func)  # TODO: figure out if these casts are sensible
        self._result_type = result_type

        requisite = {fd for fd in self._param_class.model_fields}
        kv = {k: v for k, v in params.items() if k in requisite}
        self._flagger_params = self._param_class(**kv)

    @override
    def run_process(self, run_params: RunParams) -> FlagResult:
        res = type(self)._flagger_func(run_params, self.flagger_params)
        if not isinstance(res, self._result_type):
            raise RuntimeError("Flagger returned wrong result type - engine misconfigured")
        return res

    @property
    def flagger_params(self) -> FixedParams:
        return self._flagger_params




# END SECTION: Behaviour Mixins -------------------------


class ProcessTypeEnum(Enum):
    TAGGER = ReadTaggerEngine
    FLAGGER = VariantFlaggerEngine
