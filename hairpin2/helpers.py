from enum import IntEnum, Flag


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
    print_enum: Flag
) -> None:
    print([':'.join([str(e), hex(e.value)]) for e in print_enum])


def print_enum(
    print_enum: IntEnum
) -> None:
    print([e for e in print_enum])
