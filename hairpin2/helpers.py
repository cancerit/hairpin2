# hairpin2
#
# Copyright (C) 2024 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


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


def test_options(args):
    if not (0 < args.min_clip_quality < 93):
        cleanup(msg='invalid --min-clip-quality; range 0-93')
    if not (0 < args.min_mapping_quality < 60):
        cleanup(msg='invalid --min-mapping-quality; range 0-60')
    if not (0 < args.min_base_quality < 93):
        cleanup(msg='invalid --min-base-quality; range 0-93')
    if not (0 < args.position_fraction < 1):
        cleanup(msg='invalid --position-fraction; range 0-1')


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
