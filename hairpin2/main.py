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

# pyright: reportImplicitStringConcatenation=false

from collections.abc import Generator, Iterable, Mapping
import pysam
from htsflow.rascheduler import RAExec
from hairpin2.const import VALID_NUCELOTIDES, MutTypes
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.read_preprocessors.mark_overlap import TaggerOverlap
from hairpin2.read_preprocessors.mark_support import TaggerSupporting
from htsflow.structures import FlagResult, ReadView, TestOutcomes
from hairpin2.flaggers import ADF, ALF, DVF, LQF
import logging
from itertools import tee
from typing import Any


# NOTE: requires VariantFileHeader associated with variant records to have been initialised correctly for the output
# TODO: improve on the above by yielding back the flag results not the record?
def hairpin2(
    records: Iterable[pysam.VariantRecord],  # TODO: take extended variant?
    sample_to_alignment: Mapping[str, pysam.AlignmentFile],
    configd: dict[str, Any],
    quiet: int = 0
) -> Generator[pysam.VariantRecord, Any, Any]:
    '''
    read-aware artefactual variant flagging algorithms. Flag variants in VCF using statistics calculated from supporting reads found in ALIGNMENTS, and emit the flagged VCF to stdout.
    '''

    # # test records
    proc_exec = RAExec.from_config(configd, {TaggerSupporting, TaggerOverlap, LQF.TaggerLowQual, DVF.TaggerDupmark, LQF.FlaggerLQF, ADF.FlaggerADF, ALF.FlaggerALF, DVF.FlaggerDVF})
    for record in records:
        record_flagd: dict[str, list[FlagResult]] = {}
        if record.alts is None:
            if quiet < 1:
                logging.warning(
                    f'{record.chrom: <7}:{record.pos: >12} (variant id: {record.id!r}) ¦ no alts for this record, skipping variant'
                )
        else:
            # TODO: also need to mandate/require field presence on variant records I suppose
            samples_w_mutants = [name
                                 for name
                                 in record.samples
                                 if record.samples[name]["GT"] != (0, 0)]
            if len(samples_w_mutants) == 0:
                if quiet < 1:
                    logging.warning(
                        f'{record.chrom: <7}:{record.pos: >12} (variant id: {record.id!r}) ¦ no samples exhibit alts associated with this record, skipping variant'
                    )
            else:
                reads_by_sample: dict[str, list[pysam.AlignedSegment]] = {}

                # get pileup
                for k, v in sample_to_alignment.items():
                    if k in samples_w_mutants:
                        read_iter, test_iter = tee(v.fetch(record.chrom,
                                                           record.start,
                                                           (record.start + 1)))
                        try:
                            _ = next(test_iter)
                        except StopIteration:
                            continue  # no reads for that sample cover this region  BUG: should warn!!
                        else:
                            reads_by_sample[k] = list(read_iter)

                # --- test by alt ---
                # TODO: put mutation type detection under testing
                # the ability to handle complex mutations would be a potentially interesting future feature
                # for extending to more varied artifacts
                for alt in record.alts:
                    can_parse_alt = set(alt).issubset(VALID_NUCELOTIDES)
                    alt_len = len(alt)

                    # NOTE: (intentionally) doesn't support symbolic del/ins
                    if can_parse_alt:
                        if record.rlen == alt_len:
                            mut_type = MutTypes.SUB
                        elif alt_len < record.rlen:
                            mut_type = MutTypes.DEL  # NOTE: track back how a * del would be handled
                        elif record.rlen == 1:
                            mut_type = MutTypes.INS
                        else:
                            if quiet < 1:
                                logging.warning(
                                    f'{record.chrom: <7}:{record.pos: >12} (variant id: {record.id!r}) ¦ could not infer mutation type from ALT/REF: {alt!r}/{record.ref!r}, skipping variant'
                                )
                            continue
                    else:
                        if quiet < 1:
                            logging.warning(
                                f'{record.chrom: <7}:{record.pos: >12} (variant id: {record.id!r}) ¦ ALT: {alt!r} contains nucleotides hairpin2 does not parse (only {VALID_NUCELOTIDES!r}), skipping variant'
                            )
                        continue

                    # instantiate test data obj/s
                    test_reads = ReadView(ReadView.convert_pysam_to_extread(reads_by_sample, validate=True))
                    run_data = RunParamsShared(record=record, reads=test_reads, alt=alt, mut_type=mut_type)  # TODO: allow positional args
                    for result in proc_exec.run(run_data):
                        record_flagd.setdefault(result.FlagName, []).append(result)  # if dict entry exists, append, if not, create and append

        if any(lst for lst in record_flagd.values()):
            for fname in record_flagd:
                if any(fres.variant_flagged == TestOutcomes.VARIANT_FAIL for fres in record_flagd[fname]):
                    record.filter.add(fname)
                # TODO: LQF prints no info!
                record.info.update({fname: ','.join([fl.getinfo() for fl in record_flagd[fname]])})  # pyright: ignore[reportArgumentType, reportUnknownMemberType]

        yield record

