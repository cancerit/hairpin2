# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
