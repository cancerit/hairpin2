import hairpin2.abstractflaggers as haf
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.ref2seq import ref_end_via_cigar
from pysam import AlignedSegment


OVERLAP_TAG = 'zO'


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
        overlap = check_fragment_overlap(
            read,
            run_params.record.start,
        )
        if overlap:  # if bad
            read.set_tag(OVERLAP_TAG, 1, 'i')


@haf.require_read_properties(require_tags=['MC', 'zS'])  # require support, MC tag; exclude low qual  # TODO: this ALREADY needs a make_dag type function
@haf.read_tagger(tagger_param_class=None, read_modifier_func=tag_overlap, adds_tag=OVERLAP_TAG)
class TaggerOverlap(
    haf.ReadAwareProcess  # TODO/BUG: ReadAwareProcess subclasses MUST define a specific type of run params that they use, or the contravariance with run_params is lost
): pass

