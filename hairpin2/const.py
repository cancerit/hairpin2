# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
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
from enum import StrEnum, Flag, auto
from typing import Any


class TagEnum(StrEnum):
    LOW_QUAL = "LOW-QUAL"
    SUPPORT = "SUPPORTS-VAR"
    OVERLAP = "IS-OVERLAPPING-READ2"
    STUTTER_DUP = "IS-STUTTER-DUP"


# TODO : flag names


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


DEFAULT_EXEC_CONFIG: dict[str, Any] = {
        "mark-support": {
            "enable": True,
            "require-marks": [],
            "exclude-marks": []
        },
        "mark-overlap": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": []
        },
        "mark-low-qual": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": []
        },
        "mark-duplicates": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": [
                "LOW-QUAL"
            ]
        },
        "LQF": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": [
                "IS-OVERLAPPING-READ2"
            ]
        },
        "DVF": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": [
                "IS-OVERLAPPING-READ2",
                "LOW-QUAL"
            ]
        },
        "ALF": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": [
                "LOW-QUAL",
                "IS-OVERLAPPING-READ2",
                "IS-STUTTER-DUP"
            ]
        },
        "ADF": {
            "enable": True,
            "require-marks": [
                "SUPPORTS-VAR"
            ],
            "exclude-marks": [
                "LOW-QUAL",
                "IS-OVERLAPPING-READ2",
                "IS-STUTTER-DUP"
            ]
        },
    "opts": {
        "mandate-excludes": True
    }
}

