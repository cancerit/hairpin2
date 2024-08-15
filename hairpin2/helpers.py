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


# <= - is subset of
def verify_json(jd: dict) -> bool:
    return jd.keys() <= {'vcf_in', 'vcf_out', 'alignments', 'format', 'name_mapping', 'al_filter_threshold', 'min_clip_quality', 'min_mapping_quality', 'min_base_quality', 'max_read_span', 'position_fraction'}


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
