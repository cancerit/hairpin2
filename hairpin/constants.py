from enum import Enum, IntEnum, Flag
from typing import Callable, Optional
import dataclasses as d

VERSION = '0.0.1'
EXIT_SUCCESS = 0
EXIT_FAILURE = 1

FiltCodes = IntEnum('FiltCodes',
            ['SIXTYAI', 'SIXTYBI', 'ON_THRESHOLD', 'INSUFFICIENT_READS', 'NO_MUTANTS'],
            start=0)
Ops = IntEnum('Ops',
            ['MATCH', 'INS', 'DEL', 'SKIP', 'SOFT', 'HARD', 'PAD', 'EQUAL', 'DIFF', 'BACK'],
            start = 0)
ValidatorFlags = Flag('ReadFlags',
            ['CLEAR', 'FLAG', 'MAPQUAL', 'READ_FIELDS_MISSING', 'NOT_ALIGNED', 'BAD_OP', 'NOT_ALT', 'BASEQUAL', 'SHORT', 'CLIPQUAL', 'MATE_MISSING_FIELDS', 'OVERLAP'],
            start=0)


@d.dataclass
class FilterData:
    name: str
    flag: bool = False
    code: Optional[int] = None

    def set(self):
        self.flag = True

    def __iter__(self):
        return (getattr(self, field.name) for field in d.fields(self))


@d.dataclass
class HPFilter(FilterData):
    name: str = d.field(default='HPF')

@d.dataclass
class ALFilter(FilterData):
    name: str = d.field(default='ALF')
    avg_as: Optional[float] = None

@d.dataclass
class Filters:
    AL: ALFilter
    HP: HPFilter

    def __iter__(self):
        return ((field.name, getattr(self, field.name)) for field in d.fields(self))

    def fill_field(self, field_name, value):
        if hasattr(self, field_name):
            setattr(self, field_name, value)
        else:
            raise AttributeError

    def get_field(self, field_name):
        if hasattr(self, field_name):
            return getattr(self, field_name)
        else:
            raise AttributeError


FiltReturn = Callable[..., Filters]
FlagReturn = Callable[..., int]


def print_flag(
    print_enum: Enum
) -> None:
    print([':'.join([str(e), hex(e.value)]) for e in print_enum])

def print_enum(
    print_enum: Enum
) -> None:
    print([e for e in print_enum])

