from hairpin2.const import TaggerNamespaces, Tags
from hairpin2.infrastructure.configure_funcs import make_read_processor
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import TagFragmentReads


def tag_overlap(
    run_params: RunParamsShared,
    fixed_params: None,  # pyright: ignore[reportUnusedParameter]
):
    TagFragmentReads.check_for_mates(run_params.reads.all)


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_OVERLAP,
    tagger_param_class=None,
    read_modifier_func=tag_overlap,
    adds_marks=[Tags.OVERLAP_TAG],
)
class TaggerOverlap(
    ReadAwareProcess,
):
    pass
