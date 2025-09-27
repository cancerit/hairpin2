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
from enum import Enum, StrEnum, Flag, auto
from typing import Any


class ValidNucleotides(Enum):
    """
    VCF nucleotides which hairpin2 is able to process.
    """

    ALL = set(["A", "C", "T", "G", "N", "*"])


class Strand(StrEnum):
    """
    Enum used to indicated strandedness
    """

    BOTH = "BOTH"
    F = "F"
    R = "R"


class MutTypes(StrEnum):
    """
    VCF Variant types which hairpin2 is able to process.
    """

    SUB = "SUB"
    DEL = "DEL"
    INS = "INS"


class TaggerNamespaces(StrEnum):
    """
    Process names for additive read processors (functions which add tag data to reads during execution).
    """

    MARK_SUPPORT = "mark-support"
    MARK_OVERLAP = "mark-overlap"
    MARK_LOW_QUAL = "mark-low-qual"
    MARK_STUTTER_DUP = "mark-duplicates"


class Tags(StrEnum):
    """
    Tags applied by additive read processors.
    """

    SUPPORT_TAG = "SUPPORTS-VAR"
    OVERLAP_TAG = "IS-OVERLAPPING-READ2"
    LOW_QUAL_TAG = "LOW-QUAL"
    STUTTER_DUP_TAG = "IS-STUTTER-DUP"


class FlaggerNamespaces(StrEnum):
    """
    Process/Flag names for variant flaggers, which flag variants in the output VCF.
    """

    LOW_QUAL = "LQF"
    DUPLICATION = "DVF"
    POOR_ALIGNMENT_SCORE = "ALF"
    ANOMALOUS_DISTRIBUTION = "ADF"


class ValidatorFlags(Flag):
    """
    Flags used to assess reads for support of a given variant and read quality.
    """

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
    "mark-support": {"enable": True, "require-marks": [], "exclude-marks": []},
    "mark-overlap": {"enable": True, "require-marks": ["SUPPORTS-VAR"], "exclude-marks": []},
    "mark-low-qual": {"enable": True, "require-marks": ["SUPPORTS-VAR"], "exclude-marks": []},
    "mark-duplicates": {
        "enable": True,
        "require-marks": ["SUPPORTS-VAR"],
        "exclude-marks": ["LOW-QUAL"],
    },
    "LQF": {
        "enable": True,
        "require-marks": ["SUPPORTS-VAR"],
        "exclude-marks": ["IS-OVERLAPPING-READ2"],
    },
    "DVF": {
        "enable": True,
        "require-marks": ["SUPPORTS-VAR"],
        "exclude-marks": ["IS-OVERLAPPING-READ2", "LOW-QUAL"],
    },
    "ALF": {
        "enable": True,
        "require-marks": ["SUPPORTS-VAR"],
        "exclude-marks": ["LOW-QUAL", "IS-OVERLAPPING-READ2", "IS-STUTTER-DUP"],
    },
    "ADF": {
        "enable": True,
        "require-marks": ["SUPPORTS-VAR"],
        "exclude-marks": ["LOW-QUAL", "IS-OVERLAPPING-READ2", "IS-STUTTER-DUP"],
    },
    "opts": {"mandate-excludes": True},
}
