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

For each variant examined, `haripin2` determines the mutation type (SUB, INS, DEL), fetches all reads covering the mutant position from the alignments, and then walks through the following steps.

For a more complete description of the internals (including some advanced options that hairpin2 provides), see [Advanced Usage](#advanced-usage).

#### Read Taggers
---------

These processes execute prior to flaggers, examining the reads which cover the variant position. They apply extended tags to them such that downstream processes can filter/exclude/include reads from their analyses as needed. Note that tags applied to reads are internal, and are not written back to the alignment file.

##### mark-support

`mark-support` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagSupportingReads` to tag supporting reads for the variant in question.

Reads found to support the variant are marked with the tag {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`

mark-support has no dependence on other steps

##### mark-overlap

`mark-overlap` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagFragmentOverlapReads` to tag overlapping members of fragments (read pairs).

The second read found for a given qname is marked with the tag {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`

mark-overlap only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support) - i.e. those that supporting the variant. As such, {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` is only applied to read2 of an overlapping pair which supports the variant

##### mark-low-qual
```toml
[params.mark-low-qual]
min_avg_clip_quality = 35
min_mapping_quality = 11
min_base_quality = 25
```

`mark-low-qual` uses {py:class}`~hairpin2.sci_funcs.TagLowQualReads` to mark reads that fail quality assessment.
`min_avg_clip_quality` is the minimum mean quality of aligned bases for a read which contains soft clipping.
`min_mapping_quality` is the minimum mapping quality (aka mapq) for a read.
`min_base_quality` is the minimum base quality for any position covering the variant in a read.

Reads found to be of insufficient quality are marked with the tag {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`

mark-low-qual only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. those that support the variant

##### mark-duplicates
```toml
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
---------

The following processes examine the reads covering the variant position, after they have been appropriately tagged by the read tagging processes described above. Beyond simpy flagging the variant in the FILTER field, each flagger process details the outcome of its test directly into the INFO field of the VCF record. Therefore in addition to a description of the process and relevant parameters, a description of the possible outcome reasoning is also provided. Note that any INFO field set by hairpin2 can be expanded into a more human-readable description via the [explain-var](#explain-qs-target) subtool. For all cases, the high-level outcomes are as follows:

```{list-table} **Outcomes**
:header-rows: 1
* - Name
  - Description
* - PASS
  - The variant passed the test, and therefore was not flagged.
* - FAIL
  - The variant failed the test, and therefore was flagged.
* - NA
  - The test could not be carried out.
```


##### LQF
```toml
[params.LQF]
read_loss_threshold = 0.99
min_pass_reads = 2
nsamples_threshold = 0
```

`LQF` uses {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` to flag variants which are largely supported by low quality reads.  

The test is performed by examining the proportion of reads supporting a variant marked with {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates) and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), compared to the total number of supporting reads (i.e. those tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support)).  

LQF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap) to avoid double counting.  

```{list-table} Parameters
:header-rows: 1
:widths: 30 70

* - **Name**
  - **Description**
* - `read_loss_threshold`
  - Threshold of low-quality reads against total reads. If this threshold is exceeded, there are too many low-quality reads, and the variant is flagged.
* - `min_pass_reads`
  - Minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.
```

:::{list-table} **Conditions**
:header-rows: 1

*   - Name
    - Description
*   - INSUFFICIENT_READS
    - Were there were enough reads to perform the test after read tag filtering, or not?
*   - THRESHOLD
    - Did the proportion of reads with stutter or low quality tags exceed the threshold, or not?
*   - MIN_PASS
    - Were there were enough supporting reads remaining without stutter or low quality tags, or not?
:::

In the case of a LQF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` may have `THRESHOLD`, `MIN_PASS`, or both. `NA` will only be seen with `INSUFFICIENT_READS`.


##### DVF
```toml
[params.DVF]
read_loss_threshold = 0.49
min_pass_reads = 2
nsamples_threshold = 0
```

`DVF` uses {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` to flag variants which are largely supported by reads which appear to be PCR stutter duplicates as determined by the [mark-duplicates](#mark-duplicates) process.  

The test is performed by examining the proportion of reads supporting a variant marked with {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` compared to the total number of supporting reads (i.e. those tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`). Note that this does not consider reads marked as duplicates in the alignment file itself, as it is assumed that such duplication will have been addressed by the variant caller in making the call. The aim here is to find otherwise hidden duplicates that a variant caller is unlikely to have accounted for. Also note that since the parameters are independent from [LQF](#LQF), the results of these two flags do not simply overlap.  

DVF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap) to avoid double counting, or {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual).

:::{list-table} **Parameters**
:header-rows: 1
*  - Name
   - Description
*  - read_loss_threshold
   - threshold of duplicate reads against total reads.
    If this threshold is exceeded, there are too many duplicate reads,
    and the variant is flagged.
*  - min_pass_reads
   - minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.
:::
:::{list-table} **Conditions**
:header-rows: 1

*   - Name
    - Description
*   - INSUFFICIENT_READS
    - Were there were enough reads to perform the test after read tag filtering, or not?
*   - THRESHOLD
    - Did the proportion of reads with stutter or low quality tags exceed the threshold, or not?
*   - MIN_PASS
    - Were there were enough supporting reads remaining without stutter or low quality tags, or not?
:::

In the case of a DVF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` may have `THRESHOLD`, `MIN_PASS`, or both. `NA` will only be seen with `INSUFFICIENT_READS`.

##### ALF
```toml
[params.ALF]
avg_AS_threshold = 0.93
```

`ALF` uses {py:class}`~hairpin2.sci_funcs.AlignmentScoreTest` to flag variants where the mean alignment score per base is less than `avg_AS_threshold`. The calcuation for each read is simply the value retrieved from the `AS` SAM tag divided by the read length. Users familiar with previous incarnations of LCMB variant post-processing should note that this is largely equivalent to `ASRD`.  

ALF operates only on reads tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. only reads that support the given variant.

ALF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap), or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates).

:::{list-table} **Parameters**
:header-rows: 1
*  - Name
   - Description
*  - avg_AS_threshold
   - threshold of mean alignment score per base. If this threshold is not reached, mean alignment score per base is too low, and the variant is flagged.
:::
:::{list-table} **Conditions**
:header-rows: 1

*   - Name
    - Description
*   - INSUFFICIENT_READS
    - Were there were enough reads to perform the test after read tag filtering, or not?
*   - INSUFFICIENT_AS_TAGS
    - Of the available reads for testing, do enough of them have an AS tag to perform the test, or not?
*   - ON_THRESHOLD
    - Were there were enough supporting reads remaining without stutter or low quality tags, or not?
:::

In the case of a ALF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` will only have `ON_THRESHOLD`. `NA` will be seen with one of `INSUFFICIENT_READS` or `INSUFFICIENT_AS_TAGS`.

##### ADF
```toml
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

```{admonition} From Ellis et al.:
:class: qoute

For each variant, if the number of variant-supporting reads determined is low (i.e. 0–1
reads) for one strand, follow Option 1. For each variant, if both strands have sufficient variant-supported reads
(i.e. ≥2 reads), follow Option 2.

- Option A - Low number of variant-supporting reads on one strand
    - (i) For each variant, if one strand had too few variant-supporting reads, the other strand must
    conform to:
        - Fewer than 90% of variant-supporting reads [ON THE STRAND] have the variant located within the first 15% of the read measured from the alignment start position.
        - MAD >0 [MEDIAN ABSOLUTE DEVIATION] and s.d. >4 for that strand.

- Option B - Sufficient variant-supporting reads on both strands
    - (i) For each variant, if both strands have sufficient variant-supporting reads (i.e., ≥2 reads),
    then one of the following must be true:
        - Fewer than 90% of variant-supporting reads should have the variant located within the first 15% of the read measured from the alignment start position.
        - MAD >2 and s.d. >2 for both strands.
        - MAD >1 and s.d. >10 for one strand (i.e., strong evidence of high variability in variant position in variant-supporting reads).
    
```

The point of ambiguity is whether or not, on path A, to include the single read from the strand which does not sufficiently
support the variant in the test of positional distribution across the supporting reads.
The present interpretation is that the test should be applied only to the reads on the strand which has sufficient variant-supporting reads,
since the phrasing "the other strand must conform to" implies the exclusion of the single read on the low-support
strand.  


ADF operates only on reads tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. only reads that support the given variant.  

ADF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap), or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates).

:::{list-table} **Parameters**
:header-rows: 1
*  - Name
   - Description
*  - edge_definition
   - percentage fraction of a read to be considered the edge when interpreting the conditional above, where 15% is used.  
*  - edge_clustering_threshold
   - percentage threshold for the maximum fraction of reads that may express the variant at read edge as defined by `edge_defintion` when interpreting the conditional above, where 90% is used.  
*  - min_MAD_one_strand
   - threshold for Median Absolute Deviation to be used on path A of the conditional above, where 0 is used.  
*  - min_sd_one_strand
   - threshold for standard deviation to be used on path A of the condtional above, where 2 is used.  
*  - min_MAD_both_strand_weak
   - threshold for Median Absolute Deviation to be used on path B, list option 2 of the condtional above, where 2 is used.  
*  - min_sd_both_strand_weak
   - threshold for standard deviation to be used on Path B, list option 2 of the conditional above, where 2 is used.  
*  - min_MAD_both_strand_strong
   - threshold for Median Absolute Deviation to be used on path B, list option 3 of the condtional above, where 1 is used.  
*  - min_sd_both_strand_strong
   - threshold for standard deviation to be used on Path B, list option 3 of the conditional above, where 10 is used.  
*  - low_n_supporting_reads_boundary
   - threshold for the minimum number of variant-supporting reads to be considered "low", per the first sentence of the conditional above, where 1 is used.  
*  - min_non_edge_reads
   - threshold; at least N supporting reads express variant away from the read edge as defined by edge_definition.
:::

`min_non_edge_reads` is an addition/extension from hairpin2 beyond the original test described by Ellis et al. In the default parameter configs provided, all values are set to those described in the original conditional (per the param block at the top of this section) - `min_non_edge_reads` is set to 0, effectively disabling it.

:::{list-table} **Conditions**
:header-rows: 1

*   - Name
    - Description
*   - NO_TESTABLE_READS
    - ...
*   - INSUFFICIENT_READS
    - not the same as in other flags! confusing
*   - EDGE_CLUSTERING
    - ...
*   - ONE_STRAND_DISTRIB
    - ...
*   - BOTH_STRAND_DISTRIB_BOTH
    - ...
*   - BOTH_STRAND_DISTRIB_ONE
    - ...
*   - MIN_NON_EDGE
    - ...
:::

In the case of an ADF `PASS` ...

#### Process Relationships Table

:::{list-table} **Read Taggers**
:header-rows: 1

*   - Name
    - Adds Tag
    - Requires Tags
    - Excludes Tags
*   - mark-support
    - ...
    - ...
    - ...
*   - mark-overlap
    - ...
    - ...
    - ...
*   - mark-low-qual
    - ...
    - ...
    - ...
*   - mark-duplicates
    - ...
    - ...
    - ...
:::
:::{list-table} **Variant Flaggers**
:header-rows: 1

*   - Flag
    - Requires Tags
    - Excludes Tags
*   - LQF
    - ...
    - ...
*   - DVF
    - ...
    - ...
*   - ALF
    - ...
    - ...
*   - ADF
    - ...
    - ...
:::

(understanding-decisions-target)=
### Understanding Decisions

hairpin2 stores detailed reasoning for each flagging decision in the INFO field of each tested variant. There will be VCF-spec standard key=value pairs for each flag tested. The key is the flag name, as would be applied to the FILTER column if the variant fails, e.g., `LQF`. For each flag, there will be as many key=value pairs as there are alts for the variant.  

The value side of the key=value pair is formatted as pipe (|) separated fields. The standard fields shared by all flags are as follows:

:::{list-table} **Value Fields**
*   - alt
    - The alt to which this data applies
*   - outcome
    - The result of testing for this alt (PASS, FAIL, NA)
*   - conditions
    - The conditions upon which the outcome is based, stored as a bitwise hex flag
*   - nreads
    - The number of reads examined by this test, after filtering to appropriately tagged reads
:::

A flag may add any number of subsequent extra fields. So in full, a hairpin2 INFO key=value pair looks like `<flag_name>=<alt>|<outcome>|<conditions>|<nreads>|...`.

key=value paris for any hairpin2 flag can be input to the [explain-var](#explain-qs-target) tool to convert this data into a more human-readable format, including descriptions of any additional fields. All fields are also fully described in the VCF header.

(repro-target)=
### Reproducibility & Parameter Distribution

hairpin2 stores the complete information necessary to distribute and reproduce a run directly in the VCF header. The following header keys are used

:::{list-table} **Header Keys**
*   - hairpin2_version
    - stores version of hairpin2 used
*   - hairpin2_params
    - stores paramters and execution configuration
*   - hairpin2_samples
    - stores which samples in the VCF were selected for testing
:::

The hairpin2_params key stores data in a JSON compressed with `zlib` and encoded to an ascii string using `base85`. The result is a string of ascii characters, much shorter than the original JSON, which does not invite or allow accidental editing of the parameters once written into the header (either by hand or by another tool). Since zlib compression includes an error-checking checksum, the string is guaranteed to transform back into the exact parameters encoded. Base85 encoding ensures that only VCF-safe characters are used  per the VCF format spec. Most importantly, the parameter JSON can be extracted from a VCF using [get-params](#get-params-qs-target).  


### Advanced Usage

#### Exec Config \[EXPERIMENTAL\]

The exec config serves two purposes - to allow experimentation with process interdependence, and to highlight and surface implicit process independence

