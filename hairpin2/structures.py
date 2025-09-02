# pyright: reportExplicitAny=false
# pyright: reportAny=false
from pysam import AlignedSegment
from itertools import chain
from collections.abc import Iterator, Sequence, Mapping
from typing import Any, override


class ReadView(Mapping[Any, tuple[AlignedSegment, ...]]):
    __slots__ = ("_data",)

    def __init__(self, data: Mapping[Any, Sequence[AlignedSegment]], *args: str) -> None:
        if not "no_validate" in args:
            self._validate_unique(data)  # post-init validation
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
    def _validate_unique(data: Mapping[str, Sequence[AlignedSegment]]):
        for vals in data.values():
            if len(vals) != len(set(id(el) for el in vals)):
                raise ValueError("Sequence of objects in mapping contains duplicated references pointing to the same object in memory (i.e. the same read has been included twice). Cannot create ReadView")

