# pyright: reportExplicitAny=false

## --- configuration decos


from types import new_class
from typing import Any, Callable, cast, overload
from collections.abc import Mapping, Sequence
from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.process_engines import FixedParams_T, FlagResult_T, ProcessEngineProtocol, ProcessTypeEnum, ReadTaggerEngine, RunParams_contraT, VariantFlaggerEngine


def _create_generic_configure_deco(
    process_namespace: str | None,
    factory_func: Callable[[Mapping[str, Any]], ProcessEngineProtocol],
    process_type: ProcessTypeEnum
) -> Callable[[type[ReadAwareProcess]], type[ReadAwareProcess]]:
    def deco(cls: type[ReadAwareProcess]) -> type[ReadAwareProcess]:
        if not issubclass(cls, ReadAwareProcess):  # pyright: ignore[reportUnnecessaryIsInstance]
            raise TypeError('Class does not appear to inherit from ReadAwareProcess')  # pyright: ignore[reportUnreachable]

        if cls.EngineFactory is not None:
            raise TypeError('Class already appears to be manually configured as EngineFactory is present - do not combine manual configuration with decorator')
        if cls.ProcessType is not None:
            raise TypeError('Class already appears to be manually configured as ProcessType is present - do not combine manual configuration with decorator')

        bases = (cls,)

        subclass_kwds: dict[str, Any] = {}

        if cls.ProcessNamespace is None:
            if process_namespace is None:
                raise ValueError('namespace not set in decorated class, nor provided to configuration decorator. One of these options must be used')
            subclass_ns = process_namespace
        else:
            if process_namespace is not None:
                raise ValueError('namespace set in decorated class, and provided to configuration decorator. Only one of these options may be used!')
            subclass_ns = cls.ProcessNamespace

        subclass_kwds['process_param_namespace'] = subclass_ns
        subclass_kwds['engine_factory'] = factory_func
        subclass_kwds['process_type'] = process_type

        # TODO: look at wrapt for this
        def body(ns: dict[str, Any]):
            ns["__module__"] = cls.__module__  # so user class comes from the same module
            ns["__doc__"] = getattr(cls, "__doc__", None)  # so doc is user doc (might actually want to use provided func doc)
            ns["__qualname__"] = getattr(cls, "__qualname__", cls.__name__)

        new_cls = cast(type[ReadAwareProcess], new_class(cls.__name__, bases, subclass_kwds, body))

        return new_cls
    return deco


@overload
def make_read_processor(
    *,
    tagger_param_class: None,
    read_modifier_func: Callable[[RunParams_contraT], None],
    adds_marks: Sequence[str],
    process_namespace: str | None = None
) -> Callable[[type[ReadAwareProcess]], type[ReadAwareProcess]]: ...
@overload
def make_read_processor(
    *,
    tagger_param_class: type[FixedParams_T],
    read_modifier_func: Callable[[RunParams_contraT, FixedParams_T], None],
    adds_marks: Sequence[str],
    process_namespace: str | None = None
) -> Callable[[type[ReadAwareProcess]], type[ReadAwareProcess]]: ...

def make_read_processor(
    *,
    tagger_param_class: type[FixedParams_T] | None,
    read_modifier_func: Callable[[RunParams_contraT, FixedParams_T], None] | Callable[[RunParams_contraT], None],
    adds_marks: Sequence[str],
    process_namespace: str | None = None
    # adds_read_data  # TODO
    # require_read_data  # TODO
) -> Callable[[type[ReadAwareProcess]], type[ReadAwareProcess]]:
    def init_engine(params: Mapping[str, Any]):
        return ReadTaggerEngine(
            params,
            tagger_param_class,
            read_modifier_func,
            adds_marks
        )
    return _create_generic_configure_deco(
        process_namespace,
        init_engine,
        ProcessTypeEnum.TAGGER
    )


def make_variant_flagger(
    *,
    flagger_param_class: type[FixedParams_T],
    flagger_func: Callable[[RunParams_contraT, FixedParams_T], FlagResult_T],
    result_type: type[FlagResult_T],
    process_namespace: str | None = None
    # adds_record_data ...
):
    def init_engine(params: Mapping[str, Any]):
        return VariantFlaggerEngine(
            params,
            flagger_param_class,
            flagger_func,
            result_type
        )
    return _create_generic_configure_deco(
        process_namespace,
        init_engine,
        ProcessTypeEnum.FLAGGER
    )


# NOTE: given this will only ever be created dynamically, it would be insane to use a class (as now) rather than an instance
# class DescribedFunc(ABC):
#     ProcName: ClassVar[str]
#     Func: ClassVar[Callable[..., Any]]
#     FuncSig: ClassVar[Signature]
#     _bound_args: BoundArguments | None = None

#     def __init_subclass__(
#         cls,
#         proc_name: str,
#         func: Callable[..., Any]
#     ) -> None:
#         cls.ProcName = proc_name
#         cls.Func = func
#         cls.FuncSig = inspect_sig(func)  # TODO: fail on positional only args present

#     @property
#     def requests_args(
#         self
#     ):
#         # TODO: that aren't alignedsegment/variantrecord or container of only - need to separate fixed params from runtime params
#         # ahhhh that actually might need sig explosion in the decorator then
#         # or some of these should be classmethods idk
#         return set(self.FuncSig.parameters.keys())

#     def load(
#         self,
#         params: Mapping[str, Any]
#     ):
#         kv = {k: v for k,v in params.items() if k in self.requests_args}
#         # TODO: can't bind till we have runtime params also!!
#         self._bound_args = self.FuncSig.bind(**kv)





# NOTE: actually this needs to decorate a class def since the free func needs to remain free!
# def wrap_func(func: Callable[..., Any], step_type: Literal['tagger', 'flagger'], proc_name: str):
#     sig = inspect_sig(func)


#     match step_type:
#         case 'tagger':
#             polymorphic_arg = cast(OrderedDict[str, Parameter], sig.parameters.copy()).popitem(last=False)
#             # check type is alignedsegment, list[alignedsegment] or mapping[Any, AlignedSegement], get appropriate mixin (or config) to inject based on that

#         case 'flagger':
#             # check for variant record, ReadView/Mapping/List arg - should only be one acceptable sig pattern I think. N.B. fixed args can be mutable if you need stateful info across many variants.
#             # Oh actually you could also provide a run params class to the decorator which we could look for args in first. Not important now.
#             # NOTE: trying to be too generic now will kill this project


#     bases = (DescribedFunc, )

#     subclass_kwds = {
#         'proc_name': proc_name,
#         'func': func
#     }

#     new_cls = new_class(proc_name + '_DescribedFunc', bases, subclass_kwds)
#     return new_cls
