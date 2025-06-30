import hairpin2.abstractfilters as haf
from dataclasses import field
from typing import ClassVar, override
from collections.abc import Sequence
from pysam import AlignedSegment
from enum import IntEnum, auto
from statistics import median
# pyright: reportExplicitAny=false
# pyright: reportAny=false
# # pyright: reportUnnecessaryIsInstance=false

# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


class ALCodes(IntEnum):
    INSUFFICIENT_READS = 0
    ON_THRESHOLD = auto()

class Result(haf.FilterResult[ALCodes]):
    Codes: ClassVar[type[ALCodes]] = ALCodes
    alt: str
    avg_as: float | None
    name: str = field(default='ALF', init=False)

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.avg_as}"

class Params(haf.FilterParams):
    al_thresh: float = 0.93

class Filter(haf.FilterTester[Sequence[AlignedSegment], Params, Result]):  # n.b. you can retain the broadest generic nature by reusing the TypeVar and not specifying a concrete type for any of the type arguments to FilterTester
    """
    Describe filter
    """
    @override
    def test[T: Sequence[AlignedSegment]](
        self,
        alt: str,
        variant_reads: T,
    ) -> tuple[T, Result]:
        # TODO: runtime checking of variant reads (in parent?) and params
        if len(variant_reads) < 1:
            code = Result.Codes.INSUFFICIENT_READS
            fresult = Result(
                flag=None,
                code=code,
                alt=alt,
                avg_as=None
            )
        
        aln_scores: list[float] = []

        for read in variant_reads:
            try:
                aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
            except KeyError:
                pass
        if len(aln_scores) != 0:
            avg_as = median(aln_scores)
            code = Result.Codes.ON_THRESHOLD
            flag = False
            if avg_as <= self.fixed_params.al_thresh:
                flag = True
            fresult = Result(
                flag=flag,
                code=code,
                alt=alt,
                avg_as=avg_as
            )
        else:
            code = Result.Codes.INSUFFICIENT_READS
            flag = False
            fresult = Result(
                flag=flag,
                code=code,
                alt=alt,
                avg_as=None
            )

        return variant_reads, fresult
