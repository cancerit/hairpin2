from hairpin2.infrastructure.configure_funcs import make_read_processor
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.const import Tags, TaggerNamespaces
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import ReadTaggingFuncs


def tag_overlap(
    run_params: RunParamsShared,
    fixed_params: None  # pyright: ignore[reportUnusedParameter]
):
    for read in run_params.reads.all:
        _ = ReadTaggingFuncs.check_fragment_overlap(read, run_params.record.start, mark=True)


# TODO - need some way to require and assure presence of bam-level tags, e.g. MC
@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_OVERLAP,
    tagger_param_class=None,
    read_modifier_func=tag_overlap,
    adds_marks=[Tags.OVERLAP_TAG],
)
class TaggerOverlap(
    ReadAwareProcess,  # TODO/BUG: ReadAwareProcess subclasses MUST define a specific type of run params that they use, or the contravariance with run_params is lost -- later, TODO: I don't know what I meant by this, must investigate
):
    pass
