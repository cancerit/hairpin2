from pysam import AlignedSegment
from hairpin2.abstractflaggers import ReadPreprocessor, FixedParams, RequireReadProperties
from hairpin2 import ref2seq as r2s
from typing import Any, cast, override
from hairpin2.flaggers.shared import PrefilterParamsShared, RunParamsShared
from hairpin2.readqc import qc_read


class PrefilterParamsDupmark(PrefilterParamsShared):
    pass


class FixedParamsDupmark(FixedParams):
    duplication_window_size: int = 6  # -1 disables


class RunParamsDupmark(RunParamsShared):
    pass


class DupmarkPreprocessor(
    ReadPreprocessor[PrefilterParamsDupmark, FixedParamsDupmark, RunParamsDupmark],
    RequireReadProperties,
    require_tags=['MC'],
    exclude_tags=[],
    require_fields=[]
):

    @override
    def prefilter(self):
        filtered_reads: dict[Any, list[AlignedSegment]] = {}
        for sample_key, reads in self.run_params.reads.items():
            passed_reads: list[AlignedSegment] = []
            for read in reads:
                invalid = qc_read(
                    read,
                    self.run_params.record.start,
                    self.run_params.record.stop,
                    self.run_params.alt,
                    self.run_params.mut_type,
                    self.prefilter_params.min_baseq,
                    self.prefilter_params.min_mapq,
                    self.prefilter_params.min_avg_clipq,
                    
                )
                if not invalid:
                    passed_reads.append(read)
            filtered_reads[sample_key] = passed_reads
        return filtered_reads

    @override
    def process(
        self
    ):
        if not len(self.run_params.reads.all) > 1:
            return

        for reads in self.run_params.reads.values():
            if not len(reads) > 1:
                continue

            # prep data
            sample_pair_endpoints: list[tuple[AlignedSegment, tuple[int, int, int, int]]] = []
            for read in self.run_params.reads.all:
                next_ref_end = r2s.ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
                pair_endpoints: tuple[int, int, int, int] = cast(
                tuple[int, int, int, int],
                    tuple(
                        sorted([
                            read.reference_start,
                            cast(int, read.reference_end),
                            read.next_reference_start,
                            next_ref_end
                        ])
                    )
                )
                sample_pair_endpoints.append((read, pair_endpoints))
            sample_pair_endpoints = sorted(sample_pair_endpoints, key=lambda x: x[1][0])

            # test data
            # NOTE: ideally keep highest quality read from dups, not currently done
            # NOTE: do you want to only look at supporting reads for this dupmarking? Discuss
            dup_endpoint_comparison_pool: list[tuple[int, int, int, int]] = []
            dup_endpoint_comparison_pool.append(sample_pair_endpoints[0][1])  # endpoints
            for i in range(1, len(sample_pair_endpoints)):
                read = sample_pair_endpoints[i][0]
                testing_endpoints: tuple[int, int, int, int] = sample_pair_endpoints[i][1]
                max_diff_per_comparison: list[int] = []
                for comparison_endpoints in dup_endpoint_comparison_pool:
                    endpoint_diffs: tuple[int, int, int, int] = cast(
                        tuple[int, int, int, int],
                        tuple(
                            [abs(x - y) for x, y in zip(comparison_endpoints, testing_endpoints)]
                        )
                    )
                    assert len(endpoint_diffs) == 4  # sanity check
                    max_diff_per_comparison.append(max(endpoint_diffs))
                if all([x <= self.fixed_params.duplication_window_size for x in max_diff_per_comparison]):
                    # then the read pair being examined is a duplicate of the others in the pool
                    dup_endpoint_comparison_pool.append(testing_endpoints)  # add it to the comparison pool
                    read.set_tag('zD', 1, 'i')  # TODO: surface this a bit more - this is our stutter dup tag
                    # read.is_duplicate = True  # dupmark
                else:
                    # read at i is not dup of reads in dup_endpoint_test_pool
                    # start again, test read at i against reads subsequent from i in ends_sorted
                    dup_endpoint_comparison_pool = [testing_endpoints]


