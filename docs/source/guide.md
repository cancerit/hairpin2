(guide-target)=
# Guide

```{note}
It is recommended to read [Quickstart](#quickstart-target) prior to reading the guide.
```

Given a VCF, and alignment files for the relevant samples of that VCF, `hairpin2` will return a VCF with variants flagged as follows:
- `ADF` if supporting reads express the variant with anomalous positional distribution
- `ALF` if supporting reads have low median alignment score per base
- `DVF` if largely supported by duplicate reads arising from PCR stutter.
- `LQF` if largely supported by low quality reads


### Flags

<!-- NOTE: point them to the process doc before the code -->

- `ADF`; The ADF flag is an implementation and extension of the artefact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads.  
<br>
- `ALF`; The ALF flag indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` flag – artefacts with anomalous distributions often cause a marked decrease in alignment score.  
<br>
- `DVF`; The DVF flag uses a relatively naive but effective algorithm for detecting variants that are the result of hidden duplicates, themselves a result of PCR error. In regions of low complexity, homopolymer tracts and secondary structure can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling. The DVF flag attempts to catch these variants by examining the regularity of the start and end coordinates of collections of supporting reads and their mates.  
<br>
- `LQF`; The LQF flag tests whether a read is largely supported by low quality reads (as defined by the config).  

Further details of the implementation of each flag can be found by reading the [Process and Parameters](#process-target) section below.  

All flags are tunable such that their parameters can be configured to a variety of use cases and sequencing methods.

(process-target)=
### Processes and Parameters

For an input VCF, `hairpin2` will iterate over the records therein and analyse each record via a series of interdependent scientific steps.

This section details each process and its relevant parameters (if any), in execution order. The connection/interdependence between steps is also described. The parameters are also shown in their relevant TOML table as written in a hairpin2 TOML config.

For a more complete description of the internals (including some advanced options that hairpin2 provides), see [Advanced Usage](#advanced-usage).

-----

For each variant examined, `haripin2` determines the mutation type (SUB, INS, DEL), fetches all reads covering the mutant position from the alignments, and then walks through the following steps. Note that tags applied to reads are internal, and are not written back to the alignment file.

#### Read Taggers

##### mark-support

`mark-support` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagSupportingReads` to tag supporting reads for the variant in question.

Reads found to support the variant are marked with the tag {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`

mark-support has no dependence on other steps

##### mark-overlap

`mark-overlap` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagFragmentOverlapReads` to tag overlapping members of fragments (read pairs).

The second read found for a given qname is marked with the tag {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`

mark-overlap only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support) - i.e. those that supporting the variant. As such, {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` is only applied to read2 of an overlapping pair which supports the variant

##### mark-low-qual
```
[params.mark-low-qual]
min_avg_clip_quality = 35
min_mapping_quality = 11
min_base_quality = 25
```

`mark-low-qual` uses {py:class}`~hairpin2.sci_funcs.TagLowQualReads` to mark reads that fail quality assessment.
`min_avg_clip_quality` is the minimum mean quality of aligned bases for a read which contains soft clipping.
`min_mapping_quality` is the minimum mapping quality (aka mapq) for a read.
`min_base_quality` is the minimum base quality for any position covering the variant in a read.
The test also marks reads based on their SAM flags.

Reads found to be of insufficient quality are marked with the tag {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`

mark-low-qual only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. those that support the variant

##### mark-duplicates
```
[params.mark-duplicates]
duplication_window_size = 6  # 
```

`mark-duplicates` uses {py:class}`~hairpin2.sci_funcs.TagStutterDuplicateReads` to mark reads which appear to be hidden PCR duplicates.
The algorithm is simply based on assessing the variation of endpoints amongst the tested reads. If a read start/end and it's mate start/end
are within `duplication_window_size` of the closest neighbouring read's start/end and mate start/end then those reads are considered to
be duplicates of each other.  

Reads found to be duplicates are marked with the tag {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`, excepting the highest quality read of the suspected duplicates.

mark-duplicates only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. those that support the variant.


#### Variant Flaggers

##### LQF
```
[params.LQF]
read_loss_threshold = 0.99
min_pass_reads = 2
nsamples_threshold = 0
```

`LQF` uses {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` to flag variants which are largely supported by low quality reads.  

The test is performed by examining the proportion of reads supporting a variant marked with {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates) and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), compared to the total number of supporting reads (i.e. those tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support)).  

LQF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap) to avoid double counting.  

`read_loss_threshold` is the threshold of low quality reads against total reads. If this threshold is exceeded, there are too many low quality reads, and the variant is flagged.  

`min_pass_reads` is the minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.

##### DVF
```
[params.DVF]
read_loss_threshold = 0.49
min_pass_reads = 2
nsamples_threshold = 0
```

`DVF` uses {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` to flag variants which are largely supported by reads which appear to be PCR stutter duplicates as determined by the [mark-duplicates](#mark-duplicates) process.  

The test is performed by examining the proportion of reads supporting a variant marked with {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` compared to the total number of supporting reads (i.e. those tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`). Note that this does not consider reads marked as duplicates in the alignment file itself, as it is assumed that such duplication will have been addressed by the variant caller in making the call. The aim here is to find otherwise hidden duplicates that a variant caller is unlikely to have accounted for. Also note that since the parameters are independent from [LQF](#LQF), the results of these two flags do not simply overlap.  

DVF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap) to avoid double counting, or {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual).  

`read_loss_threshold` is the threshold of duplicate reads against total reads. If this threshold is exceeded, there are too many duplicate reads, and the variant is flagged.  

`min_pass_reads` is the minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.

##### ALF
```
[params.ALF]
avg_AS_threshold = 0.93
```

`ALF` uses {py:class}`~hairpin2.sci_funcs.AlignmentScoreTest` to flag variants where the mean alignment score per base is less than `avg_AS_threshold`. The calcuation for each read is simply the value retrieved from the `AS` SAM tag divided by the read length. Users familiar with previous incarnations of LCMB variant post-processing should note that this is largely equivalent to `ASRD`.  

ALF operates only on reads tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. only reads that support the given variant.

ALF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap), or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates).

##### ADF
```
[params.ADF]
edge_definition = 0.15
edge_clustering_threshold = 0.9
min_MAD_one_strand = 0.0
min_sd_one_strand = 4.0
min_MAD_both_strand_weak = 2.0
min_sd_both_strand_weak = 2.0
min_MAD_both_strand_strong = 1.0
min_sd_both_strand_strong = 10.0
low_n_supporting_reads_boundary = 1
min_non_edge_reads = 0
```

`ADF` uses {py:class}`~hairpin2.sci_funcs.AnomalousDistributionTest` to flag variants where the variant position is anomalously distributed on supporting reads. This often indicates artefacts due to secondary structure, but can also be the result of other artefactual phenomena.

ADF is an implementation of the following conditional described by [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). The full text of the conditional is reproduced below, with editorials in []. There is a point of ambiguity in the original conditional; the interpretation that this tool has opted for is indicated by [] and is expanded upon subsequently.  

<!-- I'd love to get this into a better formatted box -->
From Ellis et al.:

---
    For each variant, if the number of variant-supporting reads determined is low (i.e. 0–1
    reads) for one strand, follow Option A. For each variant, if both strands have sufficient variant-supported reads
    (i.e. ≥2 reads), follow Option B.
    
    A) Low number of variant-supporting reads on one strand
        (i) For each variant, if one strand had too few variant-supporting reads, the other strand must
        conform to:
            - Fewer than 90% of variant-supporting reads [ON THE STRAND] have the variant located within the first 15% of the read measured from the alignment start position.
            - MAD >0 [MEDIAN ABSOLUTE DEVIATION] and s.d. >4 for that strand.
    
    B) Sufficient variant-supporting reads on both strands
        (i) For each variant, if both strands have sufficient variant-supporting reads (i.e., ≥2 reads),
        then one of the following must be true:
            - Fewer than 90% of variant-supporting reads should have the variant located within the first 15% of the read measured from the alignment start position.
            - MAD >2 and s.d. >2 for both strands.
            - MAD >1 and s.d. >10 for one strand (i.e., strong evidence of high variability in variant position in variant-supporting reads).
    
---

The point of ambiguity is whether or not, on path A, to include the single read from the strand which does not sufficiently
support the variant in the test of positional distribution across the supporting reads.
The present interpretation is that the test should be applied only to the reads on the strand which has sufficient variant-supporting reads,
since the phrasing "the other strand must conform to" implies the exclusion of the single read on the low-support
strand.  

The parameters exposed for this flag by hairpin2 are as follows:  

---

`edge_definition` - the percentage fraction of a read to be considered the edge when interpreting the conditional above, where 15% is used.  

`edge_clustering_threshold` - percentage threshold for the maximum fraction of reads that may express the variant at read edge as defined by `edge_defintion` when interpreting the conditional above, where 90% is used.  

`min_MAD_one_strand` - the threshold for Median Absolute Deviation to be used on path A of the conditional above, where 0 is used.  

`min_sd_one_strand` - the threshold for standard deviation to be used on path A of the condtional above, where 2 is used.  

`min_MAD_both_strand_weak` - the threshold for Median Absolute Deviation to be used on path B, list option 2 of the condtional above, where 2 is used.  

`min_sd_both_strand_weak` - the threshold for standard deviation to be used on Path B, list option 2 of the conditional above, where 2 is used.  

`min_MAD_both_strand_strong` - the threshold for Median Absolute Deviation to be used on path B, list option 3 of the condtional above, where 1 is used.  

`min_sd_both_strand_strong` - the threshold for standard deviation to be used on Path B, list option 3 of the conditional above, where 10 is used.  

`low_n_supporting_reads_boundary` - the threshold for the minimum number of variant-supporting reads to be considered "low", per the first sentence of the conditional above, where 1 is used.  

`min_non_edge_reads`- an additional condition is added by hairpin2 - whether at least N supporting reads express away from the read
edge. This condition is controlled by `min_non_edge_reads`.  

In the default parameter configs provided, all values are set to those described in the original conditional (per the param block at the top of this section). `min_non_edge_reads` is set to 0, effectively disabling it.

---

ADF operates only on reads tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. only reads that support the given variant.  

ADF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap), or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates).

(explain-target)=
### Understanding Decisions

explain-var TBD

(get-target)=
### Reproducibility & Parameter Distribution

param packing TBD

### Advanced Usage

#### Processes

TBD

 A process in this context is an implementation of scientific/biological logic (see {py:mod}`hairpin2.sci_funcs`) wrapped with infrastructure so `hairpin2` can expose fixed parameters via the configs, validate and filter data for each process, and so on. Each process will take data per each iteration (reads, a variant record, associated calcuations) and fixed parameters (if needed). Fixed parameters simply indicates that the value is fixed for the duration of a hairpin2 run. In addition to the infrastructure concerns described, wrapping each scientific step in these processes makes it easier to reason about the flow of execution, surface the interdependence between steps, and allows for some more ....

#### Exec Config \[EXPERIMENTAL\]

The exec config serves two purposes - to allow experimentation with process interdependence, and to highlight and surface implicit process independence

