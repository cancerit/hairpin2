from typing import Any

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
