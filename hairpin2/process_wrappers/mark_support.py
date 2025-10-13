# MIT License

# Copyright (C) 2024, 2025 Genome Research Ltd.

# Author: Alex Byrne <ab63@sanger.ac.uk>

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
