from enum import IntEnum, Flag
import logging
import sys
from hairpin2 import constants as c


def cleanup(code: int = c.EXIT_FAILURE, msg: None | str = None) -> None:
    if code != c.EXIT_SUCCESS and msg:
        logging.error(msg)
    if code == c.EXIT_SUCCESS:
        logging.info('hairpin complete')
    sys.exit(code)


def has_duplicates(
    l: list
) -> bool:
    return len(l) != len(set(l))


def lists_not_equal(
    l1: list | set,
    l2: list | set
) -> bool:
    return sorted(l1) != sorted(l2)


def print_flag(
    print_enum: Flag,
) -> None:
    pl = []
    for e in print_enum:
        vs = '-'.join([str(int(e.value)), str(hex(e.value)), str(bin(e.value))])
        pl.append(': '.join([str(e), vs]))
    print(pl)


def print_enum(
    print_enum: IntEnum
) -> None:
    print([e for e in print_enum])  # type: ignore - iterating works fine
