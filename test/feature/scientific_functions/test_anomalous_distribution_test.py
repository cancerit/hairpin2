import copy
import pytest
import pysam

from hairpin2.infrastructure.structures import TestOutcomes, ExtendedRead
from hairpin2.sci_funcs import AnomalousDistributionTest


class TestAnomalyDistributionTest:
    def test_path_insufficient_reads(self):
        """Test with no reads - should return NA due to insufficient reads"""
        anomaly_distribution_test = AnomalousDistributionTest()
        reads_in = []
        result = anomaly_distribution_test.test_variant_reads(
            reads=reads_in,
            record_start=150,
            edge_definition=0.15,
            edge_clustering_threshold=0.9,
            min_MAD_one_strand=0,
            min_sd_one_strand=4,
            min_MAD_both_strand_weak=2,
            min_sd_both_strand_weak=2,
            min_MAD_both_strand_strong=1,
            min_sd_both_strand_strong=10,
            low_n_supporting_reads_boundary=0,
            min_non_edge_reads=1,
        )
        assert result.outcome == TestOutcomes.NA
        assert result.reason == result.Info.NO_TESTABLE_READS

    def test_insufficient_reads_below_threshold(self, read: ExtendedRead):
        """Test with reads below the minimum threshold - should return NA"""
        anomaly_distribution_test = AnomalousDistributionTest()
        read_copy = copy.deepcopy(read)
        read_copy.flag = read_copy.flag | 0x10  # Reverse strand

        reads_in = [read, read_copy]
        result = anomaly_distribution_test.test_variant_reads(
            reads=reads_in,
            record_start=150,
            edge_definition=0.15,
            edge_clustering_threshold=0.9,
            min_MAD_one_strand=0,
            min_sd_one_strand=4,
            min_MAD_both_strand_weak=2,
            min_sd_both_strand_weak=2,
            min_MAD_both_strand_strong=1,
            min_sd_both_strand_strong=10,
            low_n_supporting_reads_boundary=5,
            min_non_edge_reads=1,
        )
        assert result.outcome == TestOutcomes.NA
        assert result.reason == result.Info.INSUFFICIENT_READS

