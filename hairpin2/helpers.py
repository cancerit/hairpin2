# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
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

from enum import IntEnum, Flag
import logging
import sys
from hairpin2 import constants as cnst
from typing import Any
from collections.abc import Collection
# pyright: reportExplicitAny=false


def cleanup(code: int = cnst.EXIT_FAILURE, msg: None | str = None) -> None:
    if code != cnst.EXIT_SUCCESS and msg:
        logging.error(msg)
    if code == cnst.EXIT_SUCCESS:
        logging.info('hairpin complete')
    sys.exit(code)


def has_duplicates(
    l: list[Any]
) -> bool:
    return len(l) != len(set(l))


def collection_equal(
    l1: Collection[Any],
    l2: Collection[Any]
) -> bool:
    return sorted(l1) == sorted(l2)


def print_flag(
    print_enum: Flag,
) -> None:
    pl: list[str] = []
    for e in print_enum:
        vs = '-'.join([str(int(e.value)), str(hex(e.value)), str(bin(e.value))])
        pl.append(f'{e}: {vs}')
    print(pl)


def print_enum(
    print_enum: IntEnum
) -> None:
    print([e for e in print_enum])  # pyright: ignore[reportGeneralTypeIssues] - iterating works fine
