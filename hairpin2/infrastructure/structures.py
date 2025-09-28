from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Flag, StrEnum
from pysam import AlignedSegment, VariantRecord
from itertools import chain
from collections.abc import Iterator, Sequence, Mapping
from typing import (
    ClassVar,
    Protocol,
    override,
    Any,
    Generic,
    Literal,
    TypeVar,
    final,
    overload,
    TYPE_CHECKING,
    runtime_checkable,
)


wrap_T = TypeVar("wrap_T")


class DataExtensionBase(ABC, Generic[wrap_T]):
    """
    Wrapper adding metadata handling to object, and forwarding
    other methods to the wrapped object.

    pysam objects which we expect to operate on are C extensions.
    This means no monkey patching. One could store associated metadata
    in a seperate object but the ergonomics of doing so clash with the
    intention of this library. Hence this approach - a wrapper that
    uses method forwarding to the wrapped object, and adds non-clashing
    methods of it's own to store per-object level data. Note that the
    wrapper is intentionally non-transparent, i.e. it will not pass
    an isinstance check against the type of the wrapped object.
    """

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
    def unwrapped(self):
        return self.__wrap_obj

    # NOTE: the following instance methods are all private
    # As the opinion of the package is that one should define
    # free functions specific to any extended subclass one
    # creates, encouraging greater flexibility with regard
    # to the properties of the wrapped object

    def _record_ext_op(self, op: str):
        """
        Record having performed a process on this object,
        regardless of outcome. Useful when using _ext_mark,
        as objects may be checked, but on the outcome, not
        marked. Without this record it would not be clear
        that they had been checked at all.
        """
        self.__operation_record.add(op)

    @property
    def _operation_history(self):
        return self.__operation_record

    def _ext_mark(
        self,
        mark: str,
    ):
        """
        record a truth tag on this object
        """
        self.__tag_store.add(mark)

    def _has_ext_mark(self, mark: str):
        return True if mark in self.__tag_store else False

    def _store_ext_data(self, key: str, value: Any, overwrite: bool = False):
        if key in self.__data_store and not overwrite:
            raise KeyError("Key already present and overwrite is False")
        self.__data_store[key] = value

    @overload
    def _get_ext_data[T](
        self, key: str, expected_type: type[T], raise_on_missing: Literal[True]
    ) -> T: ...
    @overload
    def _get_ext_data[T](
        self, key: str, expected_type: type[T], raise_on_missing: Literal[False]
    ) -> T | None: ...
    @overload
    def _get_ext_data(
        self, key: str, expected_type: None, raise_on_missing: Literal[True]
    ) -> Any: ...
    @overload
    def _get_ext_data(
        self, key: str, expected_type: None, raise_on_missing: Literal[False]
    ) -> Any | None: ...

    def _get_ext_data[T](
        self, key: str, expected_type: type[T] | None = None, raise_on_missing: bool = True
    ) -> T | Any | None:
        obj = self.__data_store.get(key)
        if obj is None and raise_on_missing:
            raise KeyError  # placeholder
        if expected_type is not None and not isinstance(obj, expected_type):
            raise TypeError
        return obj

    def _has_ext_data[T](
        self,
        key: str,
        # expected_type: type[T] | None = None,
    ) -> bool:
        return key in self.__data_store


def record_operation(ext_obj: DataExtensionBase[Any], op: str):
    ext_obj._record_ext_op(op)  # pyright: ignore[reportPrivateUsage]


def mark_ext_obj(ext_obj: DataExtensionBase[Any], mark: str) -> None:
    ext_obj._ext_mark(mark)  # pyright: ignore[reportPrivateUsage]


def ext_obj_has_mark(ext_obj: DataExtensionBase[Any], mark: str):
    return ext_obj._has_ext_mark(mark)  # pyright: ignore[reportPrivateUsage]


# TODO: hold private ref to Optional[ReadView] which this ExtendedRead is a part of
# record tags/operations at that level also
_AlSegShim = AlignedSegment if TYPE_CHECKING else object


@final
class ExtendedRead(  # pyright: ignore[reportUnsafeMultipleInheritance] - because we're not really inheriting
    DataExtensionBase[AlignedSegment],
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


def make_extended_read(read: AlignedSegment | ExtendedRead) -> ExtendedRead:
    """
    free func and functional style so as to minimise pain points:
        - ease of converting from AlignedSegment to ExtendedRead, no-op if already extended
        -
    """
    if isinstance(read, ExtendedRead):
        return read
    elif isinstance(read, AlignedSegment):
        return ExtendedRead(read)
    else:
        raise TypeError(
            "data extension on read objects is only supported for pysam.AlignedSegment at this time"
        )  # pyright: ignore[reportUnreachable]


# TODO: on reflection, this should attempt to mark the underlying object
# in addition to the proxy.
def mark_read(
    read: ExtendedRead,  # | AlignedSegment - TODO
    mark: str,
    # write_to_segment: bool - TODO
) -> None:
    """
    free func and functional style so as to minimise pain points:
        - since ExtendedRead uses method forwarding, type checker is unhelpful when using .method() style
        - can add behaviour on recieving unwrapped AlignedSegment rather than extended read
    """
    mark_ext_obj(read, mark)


# TODO: handle alignedsegment(?)
def read_has_mark(read: ExtendedRead, mark: str) -> bool:
    return ext_obj_has_mark(read, mark)


_VarRecShim = VariantRecord if TYPE_CHECKING else object


@final
class ExtendedVariant(
    DataExtensionBase[VariantRecord],
    _VarRecShim,
):
    pass


# TODO: record processes executed on readview object (for execution flow)
# TODO: quick access to reads by tags, both union and intersections of
Read_T = TypeVar("Read_T", ExtendedRead, AlignedSegment, covariant=True)


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
        self._data: dict[Any, tuple[Read_T, ...]] = {
            ky: tuple(vl) for ky, vl in (data or {}).items()
        }

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
    def values(self):  # pyright: ignore[reportIncompatibleMethodOverride]
        return self._data.values()

    @override
    def items(self):  # pyright: ignore[reportIncompatibleMethodOverride]
        return self._data.items()

    @property
    def all(self):
        return list(chain.from_iterable(self._data.values()))  # collate from all samples

    @staticmethod
    def validate_read_map(read_map: Mapping[str, Sequence[AlignedSegment]]):
        for vals in read_map.values():
            if len(vals) != len(set(id(el) for el in vals)):
                raise ValueError(
                    "Sequence of objects in mapping contains duplicated references pointing to the same object in memory (i.e. the same read has been included twice). Cannot create ReadView"
                )
            if not all(isinstance(rd, AlignedSegment) for rd in vals):
                raise ValueError(
                    "Sequence of objects in mapping contains objects that are not of type pysam.AlignedSegment"
                )

    @staticmethod
    def convert_pysam_to_extread(read_map: Mapping[Any, Sequence[AlignedSegment]], validate: bool):
        if validate:
            ReadView.validate_read_map(read_map)
        retd: dict[Any, list[ExtendedRead]] = {}
        for ky, vals in read_map.items():
            retd[ky] = [ExtendedRead(val) for val in vals]
        return retd


class TestOutcomes(StrEnum):
    VARIANT_PASS = "PASS"
    VARIANT_FAIL = "FAIL"
    NA = "NA"


# NOTE/TODO: FlagResult is the part of the abstractflaggers model about which I am most skeptical
# since it uses init_subclass over a decorator, making it different
# and because it requires the user to override an abstract method
@dataclass(frozen=True)
class FlagResult(ABC):
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
    InfoFlags: ClassVar[type[Flag] | None]
    InfoFlagsAllSet: ClassVar[Flag | None]
    variant_flagged: TestOutcomes
    info_flag: Flag | None

    def __post_init__(
        self,
    ):
        cls = type(self)
        if cls.InfoFlags is None and self.info_flag is not None:
            raise AttributeError(
                f"The class definition for FlagResult {cls.FlagName!r} does not specify info flags, so you cannot set an info flag"
            )
        elif cls.InfoFlags is not None and self.info_flag is None:
            raise AttributeError(
                f"The class definition for FlagResult {cls.FlagName!r} specifies info flags, so you must set an info flag"
            )
        elif self.info_flag is not None and cls.InfoFlags is not None:
            assert cls.InfoFlagsAllSet is not None
            if not self.info_flag & cls.InfoFlagsAllSet:
                raise ValueError(
                    f"Attempting to set info bits {self.info_flag}, but this FlagResult {cls.FlagName!r} only allows bits {cls.InfoFlagsAllSet}. Note that 0 is not allowed as it is reserved to indicate an unknown condition"
                )

    def __init_subclass__(
        cls, *, flag_name: str, info_enum: type[Flag] | None = None, **kwargs: Any
    ) -> None:
        if not isinstance(flag_name, str):
            raise TypeError  # pyright: ignore[reportUnreachable]
        cls.FlagName = flag_name

        if info_enum is not None:
            if not issubclass(info_enum, Flag):
                raise TypeError  # pyright: ignore[reportUnreachable]

            cls.InfoFlags = info_enum
            cls.InfoFlagsAllSet = ~info_enum(0)
        else:
            cls.InfoFlags = None
            cls.InfoFlagsAllSet = None

        # super().__init_subclass__(**kwargs)

    @abstractmethod
    def getinfo(self) -> str:
        """
        Return basic filter info in a string formatted for use in the VCF INFO field - "<flag>|<code>".

        Each filter must return INFO as it should be formatted for the VCF INFO field, or None if not applicable.
        Subclasses must override this method to return more specific info.
        """


@runtime_checkable
class ResultPackProtocol[T: Flag](Protocol):
    Info: type[T]
    outcome: TestOutcomes
    reason: T


@runtime_checkable
class VariantTestProtocol[T: ResultPackProtocol[Flag]](Protocol):
    ResultPack: type[T]

    @classmethod
    def test_variant_reads(cls, *args: Any, **kwargs: Any) -> T: ...
