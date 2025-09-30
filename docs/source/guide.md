# Guide

Given a VCF, and alignment files for the relevant samples of that VCF, `hairpin2` will return a VCF with variants flagged as follows:
- `ADF` if supporting reads express the variant with anomalous positional distribution
- `ALF` if supporting reads have low median alignment score per base
- `DVF` if largely supported by duplicate reads arising from PCR stutter.
- `LQF` if largely supported by low quality reads


### Flags

- `ADF`; The `ADF` flag is an implementation and extension of the artifact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads. See {py:class:}`~hairpin2.sci_funcs.AnomalousDistributionTest` for more info and implementation details.  
<br>
- `ALF`; The `ALF` flag indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` flag – artefacts with anomalous distributions often cause a marked decrease in alignment score. See {py:class}`~hairpin2.sci_funcs.AlignmentScoreTest`  
<br>
- `DVF`; The `DVF` flag uses a relatively naive but effective algorithm for detecting variants which are the result of PCR error. The reasoning is as follows: In regions of low complexity, homopolymer tracts and secondary structure can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling. The `DVF` flag attempts to catch these variants by examining the regularity of the start and end coordinates of collections of supporting reads and their mates. See {py:class}`~hairpin2.sci_funcs.TagStutterDuplicateReads` and {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` for more info and implementation details.  
<br>
- `LQF`; The `LQF` flag is a superset of the `DVF` flag - it tests whether a read is largely supported by both low quality reads as determined by {py:class}`~hairpin2.TagLowQualReads` and {py:class}`~hairpin2.sci_funcs.TagStutterDuplicateReads`, i.e. stutter duplicate reads are also considered low quality. Note that because the parameters for each `LQF` and `DVF` are independent, you can indepedently set the sensitivity of each - so the result of LQF is not necessarily a complete overlap with DVF (and usually is not).  


All flags are tunable such that their parameters can be configured to a variety of use cases and sequencing methods.


### Processes and Parameters

For an input VCF, `hairpin2` will iterate over the records therein and analyse each record via a series of interdependent scientific steps.

This section details each process and its relevant parameters (if any), in execution order. The connection/interdependence between steps is also described. The parameters are also shown in their relevant TOML table as written in a hairpin2 TOML config.

For a more complete description of the internals (including some advanced options that hairpin2 provides), see [Advanced Usage](#advanced usage).

-----

For each variant examined, `haripin2` determines the mutation type (SUB, INS, DEL), fetches all reads covering the mutant position from the alignments, and then walks through the following steps. Note that tags applied to reads are internal, and are not written back to the alignment file.

#### mark-support

mark-support takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagSupportingReads` to tag supporting reads for the variant in question.

Reads found to support the variant are marked with the tag "SUPPORTS-VAR"

mark-support has no dependence on other steps

#### mark-overlap

mark-overlap takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagFragmentOverlapReads` to tag overlapping members of fragments (read pairs).

The second read found for a given qname is marked with the tag "IS-OVERLAPPING-READ2"

mark-overlap only operates on reads marked with "SUPPORTS-VAR" - i.e. those that supporting the variant. As such, "IS-OVERLAPPING-READ2" is only applied to read2 of an overlapping pair which supports the variant

#### mark-low-qual
```
[params.mark-low-qual]
min_avg_clip_quality = 35
min_mapping_quality = 11
min_base_quality = 25
```

mark-low-qual uses {py:class}`~hairpin2.sci_funcs.TagLowQualReads` to mark reads that fail quality assessment.
`min_avg_clip_quality` is the minimum mean quality of aligned bases for a read which contains soft clipping.
`min_mapping_quality` is the minimum mapping quality (aka mapq) for a read.
`min_base_quality` is the minimum base quality for any position covering the variant in a read.
The test also marks reads based on their SAM flags.

Reads found to be of insufficient quality are marked with the tag "LOW-QUAL"

mark-low-qual only operates on reads marked with "SUPPORTS-VAR", i.e. those that support the variant

#### mark-duplicates
```
[params.mark-duplicates]
duplication_window_size = 6  # 
```

mark-duplicates uses {py:class}`~hairpin2.sci_funcs.TagStutterDuplicateReads` to mark reads which appear to be hidden PCR duplicates.
The algorithm is simply based on assessing the variation of endpoints amongst the tested reads. If a read start/end and it's mate start/end
are within `duplication_window_size` of the closest neighbouring read's start/end and mate start/end then those reads are considered to
be duplicates of each other.

Reads found to be duplicates are marked with the tag "IS-STUTTER-DUP", excepting the highest quality read of the suspected duplicates.

mark-duplicates only operates on reads marked with "SUPPORTS-VAR", i.e. those that support the variant.


#### LQF
```
[params.LQF]
read_loss_threshold = 0.99
min_pass_reads = 2
nsamples_threshold = 0
```

LQF uses {py:class}`~hairpin2.sci_funs.ProportionBasedTest` to flag variants which are largely supported by low quality reads.
It does so by examining the proportion of reads supporting a variant which have the "IS-STUTTER-DUP" tag and the "LOW-QUAL" tag compared to the total number of supporting reads (i.e. those tagged with "SUPPORTS-VAR")
`read_loss_threshold` is the threshold of low quality reads against total reads. If this threshold is exceeded, there are too many low quality reads, and the variant is flagged.
`min_pass_reads` is the minimum number of supporting reads without "IS-STUTTER-DUP" and "LOW-QUAL" necessary for a variant to pass the test. If fewer, a variant is flagged.

#### DVF
```
[params.DVF]
read_loss_threshold = 0.49
min_pass_reads = 2
nsamples_threshold = 0
```

#### ALF
```
[params.ALF]
avg_AS_threshold = 0.93
```

#### ADF
```
[params.ADF]
edge_definition = 0.15
edge_clustering_threshold = 0.9
min_MAD_one_strand = 0
min_sd_one_strand = 4.0
min_MAD_both_strand_weak = 2
min_sd_both_strand_weak = 2.0
min_MAD_both_strand_strong = 1
min_sd_both_strand_strong = 10.0
min_non_edge_reads = 0
low_n_supporting_reads_boundary = 1
```


### Advanced Usage

#### Processes

TBD

 A process in this context is an implementation of scientific/biological logic (see {py:mod}`hairpin2.sci_funcs`) wrapped with infrastructure so `hairpin2` can expose fixed parameters via the configs, validate and filter data for each process, and so on. Each process will take data per each iteration (reads, a variant record, associated calcuations) and fixed parameters (if needed). Fixed parameters simply indicates that the value is fixed for the duration of a hairpin2 run. In addition to the infrastructure concerns described, wrapping each scientific step in these processes makes it easier to reason about the flow of execution, surface the interdependence between steps, and allows for some more ....

#### Exec Config \[EXPERIMENTAL\]

The exec config serves two purposes - to allow experimentation with process interdependence, and to highlight and surface implicit process independence

