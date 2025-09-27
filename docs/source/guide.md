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
- `DVF`; The `DVF` flag uses a relatively naive but effective algorithm for detecting variants which are the result of PCR error. The reasoning is as follows: In regions of low complexity, homopolymer tracts and secondary structure can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling. The `DVF` flag attempts to catch these variants by examining the regularity of the start and end coordinates of collections of supporting reads and their mates. See {py:meth}`~hairpin2.sci_funcs.ReadTaggingFuncs.check_stutter_duplicates` and {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` for more info and implementation details.  
<br>
- `LQF`; The `LQF` flag is a superset of the `DVF` flag - it tests whether a read is largely supported by both low quality reads as determined by {py:meth}`~hairpin2.sci_funcs.ReadTaggingFuncs.check_low_qual_read` and {py:meth}`~hairpin2.sci_funcs.ReadTaggingFuncs.check_stutter_duplicates`, i.e. stutter duplicate reads are also considered low quality. Note that because the parameters for each `LQF` and `DVF` are independent, you can indepedently set the sensitivity of each - so the result of LQF is not necessarily a complete overlap with DVF (and usually is not).  


All flags are tunable such that their parameters can be configured to a variety of use cases and sequencing methods.


### Processes and Parameters

For an input VCF, `hairpin2` will iterate over the records therein and analyse each record via a series of interdependent scientific steps.

This section details each process and its relevant parameters (if any), in execution order. The connection/interdependence between steps is also described. The parameters are also shown in their relevant TOML table as written in a hairpin2 TOML config.

For a more complete description of the internals (including some advanced options that hairpin2 provides), see [Advanced](#advanced).

-----

#### mark-support

mark-support takes no parameters. This process uses {py:meth}`~hairpin2.sci_funcs.ReadTaggingFuncs.check_read_supporting` to tag reads as supporting the variant in question.

#### mark-overlap

mark-overlap takes no parameters. This process uses {py:meth}`~hairpin2.sci_funcs.ReadTaggingFuncs.check_fragment_overlap` to tag overlapping second members of fragments (read pairs)

#### mark-low-qual
```
[params.mark-low-qual]
min_avg_clip_quality = 35
min_mapping_quality = 11
min_base_quality = 25
```

#### mark-duplicates
```
[params.mark-duplicates]
duplication_window_size = 6
```

#### LQF
```
[params.LQF]
read_loss_threshold = 0.99
min_pass_reads = 2
nsamples_threshold = 0
```

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


### Advanced

#### Processes

TBD

 A process in this context is an implementation of scientific/biological logic (see {py:mod}`hairpin2.sci_funcs`) wrapped with infrastructure so `hairpin2` can expose fixed parameters via the configs, validate and filter data for each process, and so on. Each process will take data per each iteration (reads, a variant record, associated calcuations) and fixed parameters (if needed). Fixed parameters simply indicates that the value is fixed for the duration of a hairpin2 run. In addition to the infrastructure concerns described, wrapping each scientific step in these processes makes it easier to reason about the flow of execution, surface the interdependence between steps, and allows for some more ....
TBD


