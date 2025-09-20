from hairpin2.abstractions.configure_funcs import make_read_processor
from hairpin2.abstractions.process import ReadAwareProcess
from hairpin2.abstractions.structures import mark_read, record_operation
from hairpin2.const import TagEnum
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.utils.ref2seq import ref_end_via_cigar
from pysam import AlignedSegment


def check_fragment_overlap(
    read: AlignedSegment,
    vcf_start: int,
):
    # avoid analysing both read1 and mate if they both cover the variant
    overlap = False

    mate_cig = str(read.get_tag('MC'))  # will error if no tag

    # NOTE: introduces strand bias!!
    if read.flag & 0x80:  # if second in pair
        read_range = range(
            read.reference_start,
            read.reference_end  # pyright: ignore[reportArgumentType]
        )
        mate_range = range(
            read.next_reference_start,
            ref_end_via_cigar(
                mate_cig,
                read.next_reference_start
            )
        )
        overlapping_positions = set(read_range).intersection(mate_range)
        if vcf_start in overlapping_positions:
            overlap = True

    return overlap


def tag_overlap(
    run_params: RunParamsShared,
):
    for read in run_params.reads.all:
        if check_fragment_overlap(
            read,
            run_params.record.start,
        ):  # if overlapping second pair member
            mark_read(
                read,
                TagEnum.OVERLAP,
            )
        record_operation(read, 'mark-overlap')


# TODO - need some way to require and assure presence of bam-level tags, e.g. MC
@make_read_processor(
    process_namespace='mark-overlap',
    tagger_param_class=None,
    read_modifier_func=tag_overlap,
    adds_marks=[TagEnum.OVERLAP],
)
class TaggerOverlap(
    ReadAwareProcess,  # TODO/BUG: ReadAwareProcess subclasses MUST define a specific type of run params that they use, or the contravariance with run_params is lost -- later, TODO: I don't know what I meant by this, must investigate
): pass

