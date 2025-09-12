# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportUnnecessaryIsInstance=false
# from abc import ABC
from abc import ABC
from pysam import AlignedSegment, VariantRecord
from itertools import chain
from collections.abc import Iterator, Sequence, Mapping
from typing import override, Any, ClassVar, Generic, Literal, TypeVar, TypedDict, final, overload, TYPE_CHECKING 



# Ideally I'd like to pass a typeddict or similar with a fixed set of keys
# such that get ext data has stronger return guarantees
wrap_T = TypeVar("wrap_T")
class WrappedWithDataStore(ABC, Generic[wrap_T]):
    __wrap_obj: wrap_T
    __data_store: dict[str, Any]
    __tag_store: set[str]
    __operation_record: set[str]

    def __init__(
        self,
        wrap_obj: wrap_T,  # I guess I trust you
    ):
        self.__wrap_obj = wrap_obj
        self.__data_store = {}
        self.__tag_store = set()
        self.__operation_record = set()

    # pass through to wrapped
    def __getattr__(self, name: str) -> Any:
        return getattr(self.__wrap_obj, name)

    @property
    def unwrapped(
        self
    ):
        return self.__wrap_obj

    def record_ext_op(
        self,
        op: str
    ):
        self.__operation_record.add(op)

    @property
    def get_operation_history(
        self
    ):
        return self.__operation_record

    def ext_mark(
        self,
        mark: str,
    ):
        self.__tag_store.add(mark)

    def has_ext_mark(
        self,
        mark: str
    ):
        return True if mark in self.__tag_store else False

    def store_ext_data(
        self,
        key: str,
        value: Any,
        overwrite: bool = False
    ):
        if key in self.__data_store and not overwrite:
            raise KeyError('Key already present and overwrite is False')
        self.__data_store[key] = value

    @overload
    def get_ext_data[T](
        self,
        key: str,
        expected_type: type[T],
        raise_on_missing: Literal[True]
    ) -> T: ...
    @overload
    def get_ext_data[T](
        self,
        key: str,
        expected_type: type[T],
        raise_on_missing: Literal[False]
    ) -> T | None: ...
    @overload
    def get_ext_data(
        self,
        key: str,
        expected_type: None,
        raise_on_missing: Literal[True]
    ) -> Any: ...
    @overload
    def get_ext_data(
        self,
        key: str,
        expected_type: None,
        raise_on_missing: Literal[False]
    ) -> Any | None: ...

    def get_ext_data[T](
        self,
        key: str,
        expected_type: type[T] | None = None,
        raise_on_missing: bool = True
    ) -> T | Any | None:
        obj = self.__data_store.get(key)
        if obj is None and raise_on_missing:
            raise KeyError  # placeholder
        if expected_type is not None and not isinstance(obj, expected_type):
            raise TypeError
        return obj


_AlSegShim = AlignedSegment if TYPE_CHECKING else object
@final
class ExtendedRead(  # pyright: ignore[reportUnsafeMultipleInheritance] - because we're not really inheriting
    WrappedWithDataStore[AlignedSegment],
    _AlSegShim,
):

    # TODO
    # def write_mark_to_tag(
    #     self,
    #     key: str,
    #     val: Any | None = None,
    #     sam_tag_type: str | None,
    #     tag_alias: str | None = None
    # ):
    #     val = self.has_ext_mark(key)
    #     assert isinstance(val, bool)
    #     sam_alias = key
    #     if tag_alias is None:
    #         if len(sam_alias) > 2:
    #             raise ValueError('sam tags may only be two characters long - use tag_alias')
    #     if tag_alias:
    #         if len(tag_alias) > 2:
    #             raise ValueError('sam tags may only be two characters long')
    #         sam_alias = tag_alias
    #     self.unwrapped.set_tag(sam_alias, val, 'i')
    pass


_VarRecShim = VariantRecord if TYPE_CHECKING else object
@final
class ExtendedVariant(
    WrappedWithDataStore[VariantRecord],
    _VarRecShim,
): pass


# TODO: record processes executed on readview object (for execution flow)
# TODO: quick access to reads by tags, both union and intersections of
# TODO: Generic over ExtendedRead, AlignedSegment

Read_T = TypeVar('Read_T', ExtendedRead, AlignedSegment, covariant=True)
class ReadView(Mapping[Any, tuple[Read_T, ...]]):
    __slots__ = ("_data",)  # pyright: ignore[reportUnannotatedClassAttribute]

    def __init__(
        self,
        data: Mapping[Any, Sequence[Read_T]],
        # *,
        # _internal_switches: list[str] | None = None  # hidden dev options
    ) -> None:
        # switches = _internal_switches or []
        # create shallow private copy with frozen (tuple) sequences
        self._data: dict[Any, tuple[Read_T, ...]] = {ky: tuple(vl) for ky, vl in (data or {}).items()}

    # mapping interface with no modification surface
    @override
    def __getitem__(self, key: Any) -> tuple[Read_T, ...]:
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
    def validate_read_map(read_map: Mapping[str, Sequence[AlignedSegment]]):
        for vals in read_map.values():
            if len(vals) != len(set(id(el) for el in vals)):
                raise ValueError("Sequence of objects in mapping contains duplicated references pointing to the same object in memory (i.e. the same read has been included twice). Cannot create ReadView")
            if not all(isinstance(rd, AlignedSegment) for rd in vals):
                raise ValueError("Sequence of objects in mapping contains objects that are not of type pysam.AlignedSegment")

    @staticmethod
    def convert_pysam_to_extread(
        read_map: Mapping[Any, Sequence[AlignedSegment]],
        validate: bool
    ):
        if validate:
            ReadView.validate_read_map(read_map)
        retd: dict[Any, list[ExtendedRead]] = {}
        for ky, vals in read_map.items():
            retd[ky] = [ExtendedRead(val) for val in vals]
        return retd

