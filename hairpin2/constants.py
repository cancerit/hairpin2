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
from typing import Callable
import dataclasses as d

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

DEFAULTS: dict[str, int | float] = dict((('al_filter_threshold', 0.93),
                                         ('min_clip_quality', 35),
                                         ('min_mapping_quality', 11),
                                         ('min_base_quality', 25),
                                         ('max_read_span', 6),
                                         ('edge_definition', 0.15),
                                         ('edge_fraction', 0.9),
                                         ('min_MAD_one_strand', 0),
                                         ('min_sd_one_strand', 4),
                                         ('min_MAD_both_strand_weak', 2),
                                         ('min_sd_both_strand_weak', 2),
                                         ('min_MAD_both_strand_strong', 1),
                                         ('min_sd_both_strand_strong', 10),
                                         ('min_reads', 1)))

FiltCodes = IntEnum('FiltCodes',
                    ['SIXTYAI',
                     'SIXTYBI',
                     'ON_THRESHOLD',
                     'INSUFFICIENT_READS',
                     'NO_MUTANTS'],
                    start=0)
Ops = IntEnum('Ops',
              ['MATCH',
               'INS',
               'DEL',
               'SKIP',
               'SOFT',
               'HARD',
               'PAD',
               'EQUAL',
               'DIFF',
               'BACK'],
              start=0)
ValidatorFlags = Flag('ValidatorFlags',
                      ['CLEAR',
                       'FLAG',
                       'MAPQUAL',
                       'READ_FIELDS_MISSING',
                       'NOT_ALIGNED',
                       'BAD_OP',
                       'NOT_ALT',
                       'BASEQUAL',
                       'SHORT',
                       'CLIPQUAL',
                       'OVERLAP'],
                      start=0)


class NoAlts(ValueError):
    pass


class NoMutants(ValueError):
    pass


@d.dataclass
class FilterData:
    name: str
    flag: bool = False
    code: int | None = None

    def set(self):
        self.flag = True

    def __iter__(self):
        return (getattr(self, field.name) for field in d.fields(self))


@d.dataclass
class ADFilter(FilterData):
    name: str = d.field(default='ADF')


@d.dataclass
class ALFilter(FilterData):
    name: str = d.field(default='ALF')
    avg_as: float | None = None


@d.dataclass
class Filters:
    AL: ALFilter
    HP: ADFilter

    def __iter__(self):
        return (getattr(self, field.name) for field in d.fields(self))

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
