import hairpin2.abstractfilters as haf
from pydantic import Field
from pydantic.dataclasses import dataclass
from typing import override
from collections.abc import Sequence
from pysam import AlignedSegment
from enum import IntEnum, auto
from statistics import median
# If you're here just to examine the scientific implementation of each filter,
# examine the `test` methods for each one
# the rest is largely boilerplate/typing magic to make the filter implementation modular and robust


class ALCodes(IntEnum):
    INSUFFICIENT_READS = 0
    INSUFFICIENT_AS_TAGS = auto()
    ON_THRESHOLD = auto()


class Result(haf.FilterResult[ALCodes]):
    alt: str
    avg_as: float | None
    name: str = Field(default='ALF', init=False)

    @override
    def getinfo(self) -> str:
        return f"{self.alt}|{self.flag}|{self.code}|{self.avg_as}"


@dataclass(frozen=True, slots=True)
class Params(haf.FilterParams):
    al_thresh: float = 0.93


class Filter(haf.FilterTester[Sequence[AlignedSegment], Params, Result]):
    """
    Describe filter
    """
    @override
    def test[T: Sequence[AlignedSegment]](
        self,
        alt: str,
        variant_reads: T,
    ) -> tuple[T, Result]:
        if len(variant_reads) < 1:
            code = ALCodes.INSUFFICIENT_READS
            fresult = Result(
                flag=None,
                code=code,
                alt=alt,
                avg_as=None
            )
        else:
            aln_scores: list[float] = []

            for read in variant_reads:
                try:
                    aln_scores.append(int(read.get_tag('AS')) / read.query_length)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]  TODO: look into fixing pysam typing
                except KeyError:
                    pass
            if len(aln_scores) != 0:
                avg_as = median(aln_scores)
                code = ALCodes.ON_THRESHOLD
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
                code = ALCodes.INSUFFICIENT_AS_TAGS
                flag = None
                fresult = Result(
                    flag=flag,
                    code=code,
                    alt=alt,
                    avg_as=None
                )

        return variant_reads, fresult
