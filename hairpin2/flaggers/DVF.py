# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import hairpin2.abstractflaggers as haf
from pydantic.dataclasses import dataclass
from typing import ClassVar, override
from enum import IntEnum, auto
from hairpin2.flaggers.shared import AltVarParams


class CodesDVF(IntEnum):
    INSUFFICIENT_READS = 0
    DUPLICATION = auto()


class ResultDVF(haf.FilterResult[CodesDVF]):
    Name: ClassVar[str] = 'DVF'
    alt: str
    loss_ratio: float  # 0 == no loss

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.loss_ratio}"  # TODO: report which samples?


class VarParamsDVF(AltVarParams): pass


@dataclass(slots=True, frozen=True)
class FixedParamsDVF(haf.FixedParams):
    duplication_window_size: int = 6  # -1 disables
    # NOTE: neither of these options prevent read removal due to duplication,
    # so `DVF.test()` always functions as QC and may still drop reads
    # - I think this is fine, just document more
    read_loss_threshold: float = 0.49  # percent threshold of N duplicate reads compared to N input reads for a given variant and sample, above which we call DVF
    nsamples_threshold: int = 0  # TODO: I'm not sure this param makes sense. I guess in a multi sample VCF it would imply less confidence in the call if only 1 sample reported duplication. But you'd still probably want to know about that sample? Discuss with Peter


class FlaggerDVF(haf.Flagger[FixedParamsDVF, VarParamsDVF, ResultDVF]):
    """
    duplication variant filter - a portion of the reads supporting the variant
    are suspected to arise from duplicated reads that have escaped dupmarking.

    In regions of low complexity, short repeats and homopolymer tracts can cause PCR stuttering.
    Leading to, for example, an additional A on the read when amplifying a tract of As.
    If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference
    between reads that are in fact duplicates. Because of this, these duplicates then evade dupmarking and give rise to
    spurious variants when calling.

    `min_boundary_deviation` sets the minimum deviation start/end coordinates, above which reads are assumed not to be duplicated
    `read_number_difference_threshold` sets the the threshold for absolute difference between the number of reads supporting the variant
    with and without duplicates removed. If this threshold is exceeded, the flag will be set.
    """
    # detect PCR duplicates previously missed due to slippage
    # this implementation assumes that sorting on first element of each sublist
    # is appropriate, per Peter's initial implementation.
    # is an all against all comparison between all read lists more appropriate?
    # and between pairs of readlists, why is comparing sorted pairs most appropriate?
    @override
    def test(
        self,
    ) -> ResultDVF:
        """
        A naive algorithm using start/end co-ordinates of read pairs to identify likely stutter duplicate reads missed by traditional dupmarking.
        """
        # NOTE: surely this shouldn't be across samples... I don't know, maybe?
        nsamples_with_duplication = 0
        nreads_by_sample: dict[str, int] = { k: len(v) for k, v in self.var_params.reads.items() }
        loss_ratio: list[float] = []

        if not any([nreads > 1 for nreads in nreads_by_sample.values()]):
            fresult = ResultDVF(
                flag=None,
                code=CodesDVF.INSUFFICIENT_READS,
                alt=self.var_params.alt,
                loss_ratio=0
            )
        else:
            code = CodesDVF.DUPLICATION  # testing possible, and this is the only relevant code
            for sample_key, reads in self.var_params.reads.items():
                ntotal = nreads_by_sample[sample_key]
                sample_loss_ratio = 0
                if ntotal > 1:
                    ndup = sum(read.is_duplicate for read in reads)
                    ntrue = abs(ndup - ntotal)
                    assert ntotal > ndup > -1  # sanity check - TODO: should probably throw an interpretable error
                    assert ntotal >= ntrue > -1
                    sample_loss_ratio = ndup / ntotal
                    if sample_loss_ratio > self.fixed_params.read_loss_threshold or ntrue < 2:
                        nsamples_with_duplication += 1

                loss_ratio.append(sample_loss_ratio)

            if nsamples_with_duplication > self.fixed_params.nsamples_threshold:
                flag = True
            else:
                flag = False
            fresult = ResultDVF(
                flag=flag,
                code=code,
                alt=self.var_params.alt,
                loss_ratio=sum(loss_ratio) / len(loss_ratio)  # TODO: discuss whether averaging is the best choice
            )

        return fresult

