# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportUnnecessaryIsInstance=false
# from abc import ABC
from pysam import AlignedSegment  #, VariantRecord
from itertools import chain
from collections.abc import Iterator, Sequence, Mapping
from typing import override, Any  #, ClassVar, Generic, Literal, TypeVar, TypedDict, final, overload, TYPE_CHECKING 



# Ideally I'd like to pass a typeddict or similar with a fixed set of keys
# such that get ext data has stronger return guarantees
# wrap_T = TypeVar("wrap_T")
# class WrappedWithDataStore(ABC, Generic[wrap_T]):
#     __wrap_obj: wrap_T
#     __data_store: dict[str, Any]

#     def __init__(
#         self,
#         wrap_obj: wrap_T,  # I guess I trust you
#         **kwargs: Any
#     ):
#         self.__wrap_obj = wrap_obj
#         self.__data_store = {}

#     # pass through to wrapped
#     def __getattr__(self, name: str) -> Any:
#         return getattr(self.__wrap_obj, name)

#     def store_ext_data(
#         self,
#         key: str,
#         value: Any
#     ):
#         self.__data_store[key] = value

#     @overload
#     def get_ext_data[T](
#         self,
#         key: str,
#         expected_type: type[T],
#         raise_on_missing: Literal[True]
#     ) -> T: ...
#     @overload
#     def get_ext_data[T](
#         self,
#         key: str,
#         expected_type: type[T],
#         raise_on_missing: Literal[False]
#     ) -> T | None: ...
#     @overload
#     def get_ext_data(
#         self,
#         key: str,
#         expected_type: None,
#         raise_on_missing: Literal[True]
#     ) -> Any: ...
#     @overload
#     def get_ext_data(
#         self,
#         key: str,
#         expected_type: None,
#         raise_on_missing: Literal[False]
#     ) -> Any | None: ...

#     def get_ext_data[T](
#         self,
#         key: str,
#         expected_type: type[T] | None = None,
#         raise_on_missing: bool = True
#     ) -> T | Any | None:
#         obj = self.__data_store.get(key)
#         if obj is None and raise_on_missing:
#             raise KeyError  # placeholder
#         if expected_type is not None and not isinstance(obj, expected_type):
#             raise TypeError
#         return obj


# _AlSegShim = AlignedSegment if TYPE_CHECKING else object
# @final
# class DataStoreRead(  # pyright: ignore[reportUnsafeMultipleInheritance] - because we're not really inheriting
#     WrappedWithDataStore[AlignedSegment],
#     _AlSegShim,
# ): pass


# _VarRecShim = VariantRecord if TYPE_CHECKING else object
# @final
# class DataStoreVariant(
#     WrappedWithDataStore[VariantRecord],
#     _VarRecShim,
# ): pass


# TODO: record processes executed on readview object (for execution flow)
# TODO: quick access to reads by tags, both union and intersections of
class ReadView(Mapping[Any, tuple[AlignedSegment, ...]]):
    __slots__ = ("_data",)  # pyright: ignore[reportUnannotatedClassAttribute]

    def __init__(
        self,
        data: Mapping[Any, Sequence[AlignedSegment]],
        *,
        _internal_switches: list[str] | None = None  # hidden dev options
    ) -> None:
        switches = _internal_switches or []
        if "no_validate" not in switches:
            self._validate(data)  # post-init validation
        # create shallow private copy with frozen (tuple) sequences
        self._data: dict[Any, tuple[AlignedSegment, ...]] = {ky: tuple(vl) for ky, vl in (data or {}).items()}

    # mapping interface with no modification surface
    @override
    def __getitem__(self, key: Any):
        return self._data[key]

    @override
    def __iter__(self) -> Iterator[Any]:
        return iter(self._data)

    @override
    def __len__(self) -> int:
        return len(self._data)

    @override
    def keys(self):
        return self._data.keys()

    @override
    def values(self):
        return self._data.values()

    @override
    def items(self):
        return self._data.items()

    @property
    def all(
        self
    ):
        return list(chain.from_iterable(self._data.values()))  # collate from all samples

    @staticmethod
    def _validate(data: Mapping[str, Sequence[AlignedSegment]]):
        for vals in data.values():
            if len(vals) != len(set(id(el) for el in vals)):
                raise ValueError("Sequence of objects in mapping contains duplicated references pointing to the same object in memory (i.e. the same read has been included twice). Cannot create ReadView")
            if not all(isinstance(rd, AlignedSegment) for rd in vals):
                raise ValueError("Sequence of objects in mapping contains objects that are not of type pysam.AlignedSegment")

