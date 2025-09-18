# pyright: reportExplicitAny=false
# pyright: reportAny=false
# pyright: reportUnnecessaryIsInstance=false
# from abc import ABC
from abc import ABC, abstractmethod
from pydantic import BaseModel, ConfigDict, field_validator
from pysam import AlignedSegment, VariantRecord
from itertools import chain
from collections.abc import Iterable, Iterator, Sequence, Mapping
from typing import ClassVar, cast, override, Any, Generic, Literal, TypeVar, final, overload, TYPE_CHECKING 


# Ideally I'd like to pass a typeddict or similar with a fixed set of keys
# such that get ext data has stronger return guarantees
wrap_T = TypeVar("wrap_T")
class _DataExtensionBase(ABC, Generic[wrap_T]):
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

    def _record_ext_op(
        self,
        op: str
    ):
        self.__operation_record.add(op)

    @property
    def _operation_history(
        self
    ):
        return self.__operation_record

    def _ext_mark(
        self,
        mark: str,
    ):
        self.__tag_store.add(mark)

    def _has_ext_mark(
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
    _DataExtensionBase[AlignedSegment],
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


def make_extended_read(
    read: AlignedSegment | ExtendedRead
) -> ExtendedRead:
    if isinstance(read, ExtendedRead):
        return read
    elif isinstance(read, AlignedSegment):
        return ExtendedRead(read)
    else:
        raise TypeError('data extension on read objects is only supported for pysam.AlignedSegment at this time')  # pyright: ignore[reportUnreachable]


def mark_read(
    read: ExtendedRead,
    mark: str
) -> None:
    read._ext_mark(mark)  # pyright: ignore[reportPrivateUsage]


_VarRecShim = VariantRecord if TYPE_CHECKING else object
@final
class ExtendedVariant(
    _DataExtensionBase[VariantRecord],
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
    def values(self):  # pyright: ignore[reportIncompatibleMethodOverride]
        return self._data.values()

    @override
    def items(self):  # pyright: ignore[reportIncompatibleMethodOverride]
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

    model_config: ConfigDict = ConfigDict(  # pyright: ignore[reportIncompatibleVariableOverride]
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
        result_codes: Iterable[int] | Sequence[str] | None,
        **kwargs: Any
    ) -> None:

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


