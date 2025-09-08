import hairpin2.abstractions.readawareproc as haf
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.const import ValidatorFlags
from hairpin2.utils.ref2seq import ref2querypos
from pysam import AlignedSegment
from typing import Literal, cast


SUPPORT_TAG = 'zS'


def check_read_supporting(
    read: AlignedSegment,
    mut_type: Literal['S', 'I', 'D'],
    alt: str,
    vcf_start: int,
    vcf_stop: int,
) -> ValidatorFlags:
    if mut_type not in ['S', 'D', 'I']:
        raise ValueError(
            'unsupported mut_type: {} - supports \'S\' (SUB) \'D\' (DEL) \'I\' (INS)'.format(mut_type))

    invalid_flag = ValidatorFlags.CLEAR
    if mut_type in ['S', 'I']:
        try:
            mut_pos = ref2querypos(read, vcf_start)
        except ValueError:
            invalid_flag |= ValidatorFlags.NOT_ALIGNED
        else:
            if mut_type == 'S':  # SUB
                if cast(str, read.query_sequence)[mut_pos:mut_pos + len(alt)] != alt:
                    invalid_flag |= ValidatorFlags.NOT_ALT
            if mut_type == 'I':  # INS - mut_pos is position immediately before insertion
                if mut_pos + len(alt) > read.query_length:
                    invalid_flag |= ValidatorFlags.SHORT
                else:
                    mut_alns = [(q, r)
                                for q, r
                                in read.get_aligned_pairs()
                                if q in range(mut_pos + 1, mut_pos + len(alt) + 1)]
                    if any([r is not None for _, r in mut_alns]):
                        invalid_flag |= ValidatorFlags.BAD_OP
                    if cast(str, read.query_sequence)[mut_pos + 1:mut_pos + len(alt) + 1] != alt:
                        invalid_flag |= ValidatorFlags.NOT_ALT
    elif mut_type == 'D':  # DEL
        rng = list(range(vcf_start, vcf_stop + 1))
        mut_alns = [q
                    for q, r
                    in read.get_aligned_pairs()
                    if r in rng]
        if len(mut_alns) != len(rng):
            invalid_flag |= ValidatorFlags.SHORT
        if (any([x is not None for x in mut_alns[1:-1]]) or
            any([x is None for x in [mut_alns[0], mut_alns[-1]]])):
                invalid_flag |= ValidatorFlags.BAD_OP

    return invalid_flag


def tag_supporting(
    run_params: RunParamsShared,
):
    for read in run_params.reads.all:
        if check_read_supporting(
            read,
            run_params.mut_type,
            run_params.alt,
            run_params.record.start,
            run_params.record.stop,
        ) == ValidatorFlags.CLEAR:  # if good
            read.set_tag(SUPPORT_TAG, 1, 'i')


@haf.read_tagger(tagger_param_class=None, read_modifier_func=tag_supporting, adds_tag=SUPPORT_TAG)  # TODO: guard against requiring/excluding and adding the same tag!
class TaggerSupporting(
    haf.ReadAwareProcess,  # TODO/BUG: ReadAwareProcess subclasses MUST define a specific type of run params that they use, or the contravariance with run_params is lost
    process_name='mark-support'
): pass

