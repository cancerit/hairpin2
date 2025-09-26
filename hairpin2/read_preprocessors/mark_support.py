from htsflow.configure_funcs import make_read_processor
from htsflow.process import ReadAwareProcess
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.const import Tags, TaggerNamespaces
from hairpin2.sci_funcs import check_read_supporting


def tag_supporting(
    run_params: RunParamsShared,
):
    for read in run_params.reads.all:
        _ = check_read_supporting(
            read,
            run_params.mut_type,
            run_params.alt,
            run_params.record.start,
            run_params.record.stop,
            mark=True
        )


@make_read_processor(
    process_namespace=TaggerNamespaces.MARK_SUPPORT,
    tagger_param_class=None,
    read_modifier_func=tag_supporting,
    adds_marks=[Tags.SUPPORT_TAG],
)
class TaggerSupporting(
    ReadAwareProcess,
): pass

