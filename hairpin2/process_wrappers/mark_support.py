from hairpin2.const import TaggerNamespaces, Tags
from hairpin2.infrastructure.configure_funcs import make_read_processor
from hairpin2.infrastructure.process import ReadAwareProcess
from hairpin2.process_wrappers.shared import RunParamsShared
from hairpin2.sci_funcs import TagSupportingReads


def tag_supporting(
    run_params: RunParamsShared,
    fixed_params: None,  # pyright: ignore[reportUnusedParameter]
):
    for read in run_params.reads.all:
        _ = TagSupportingReads.check_read_supporting(
            read,
            run_params.mut_type,
            run_params.alt,
            run_params.record.start,
            run_params.record.stop,
            mark=True,
        )


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_SUPPORT,
    tagger_param_class=None,
    read_modifier_func=tag_supporting,
    adds_marks=[Tags.SUPPORT_TAG],
)
class TaggerSupporting(
    ReadAwareProcess,
):
    pass
