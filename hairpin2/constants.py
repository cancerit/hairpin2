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
from enum import IntEnum, Flag, auto

EXIT_SUCCESS = 0
EXIT_FAILURE = 1

DEFAULTS: dict[str, int | float] = dict((('al_filter_threshold', 0.93),
                                         ('min_clip_quality', 35),
                                         ('min_mapping_quality', 11),
                                         ('min_base_quality', 25),
                                         ('duplication_window_size', 6),
                                         ('edge_definition', 0.15),
                                         ('edge_fraction', 0.9),
                                         ('min_MAD_one_strand', 0),
                                         ('min_sd_one_strand', 4),
                                         ('min_MAD_both_strand_weak', 2),
                                         ('min_sd_both_strand_weak', 2),
                                         ('min_MAD_both_strand_strong', 1),
                                         ('min_sd_both_strand_strong', 10),
                                         ('min_reads', 1)))


class Ops(IntEnum):
    MATCH = 0
    INS = auto()
    DEL = auto()
    SKIP = auto()
    SOFT = auto()
    HARD = auto()
    PAD = auto()
    EQUAL = auto()
    DIFF = auto()
    BACK = auto()

class ValidatorFlags(Flag):
    CLEAR = 0
    FLAG = auto()
    MAPQUAL = auto()
    READ_FIELDS_MISSING = auto()
    NOT_ALIGNED = auto()
    BAD_OP = auto()
    NOT_ALT = auto()
    BASEQUAL = auto()
    SHORT = auto()
    CLIPQUAL = auto()
    OVERLAP = auto()


class NoAlts(ValueError):
    pass


class NoMutants(ValueError):
    pass

