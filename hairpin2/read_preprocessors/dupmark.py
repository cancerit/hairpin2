from pysam import AlignedSegment
from collections.abc import Sequence
from hairpin2 import ref2seq as r2s
from typing import cast

def mark_stutter_dups(
    reads: Sequence[AlignedSegment],  # modifies reads in place
    dup_window_size: int
):
    if not len(reads) > 1:
        return

    if len(reads) != len(set(id(rd) for rd in reads)):
        raise ValueError("Sequence of objects submitted to function contains duplicated references pointing to the same object in memory (i.e. the same read has been included twice)")

    # prep data
    sample_pair_endpoints: list[tuple[AlignedSegment, tuple[int, int, int, int]]] = []
    for read in reads:
        if read.reference_end is None:
            raise ValueError(f"read {read.query_name} has no reference end coordinate")
        next_ref_end = r2s.ref_end_via_cigar(str(read.get_tag('MC')), read.next_reference_start)  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        pair_endpoints: tuple[int, int, int, int] = cast(
        tuple[int, int, int, int],
            tuple(
                sorted([
                    read.reference_start,
                    read.reference_end,
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
        if all([x <= dup_window_size for x in max_diff_per_comparison]):
            # then the read pair being examined is a duplicate of the others in the pool
            dup_endpoint_comparison_pool.append(testing_endpoints)  # add it to the comparison pool
            read.is_duplicate = True  # dupmark
        else:
            # read at i is not dup of reads in dup_endpoint_test_pool
            # start again, test read at i against reads subsequent from i in ends_sorted
            dup_endpoint_comparison_pool = [testing_endpoints]
