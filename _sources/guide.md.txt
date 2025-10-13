(guide-target)=
# Guide

```{note}
It is recommended to read [Quickstart](#quickstart-target) prior to reading the guide.
```

Given a VCF, and alignment files for the relevant samples of that VCF, `hairpin2` will return a VCF with variants flagged in the FILTER column as follows:

:::{list-table}
:width: 100
:widths: 15 25 60
:header-rows: 1

*   - FILTER
    - Expanded Name
    - Description
*   - **ADF**
    - **A**nomalous **D**istribution **F**lag
    - An implementation and extension of the artefact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads.
*   - **ALF**
    - **AL**ignment score **F**lag
    - Indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` flag – artefacts with anomalous distributions often cause a marked decrease in alignment score.
*   - **DVF**
    - **D**uplication **V**ariant **F**lag
    -  Uses a relatively naive but effective algorithm for detecting variants that are the result of [hidden duplicates](#dup-explain-target), themselves a result of PCR error.
*   - **LQF**
    - **L**ow **Q**uality **F**lag
    - Tests whether a variant is largely supported by low quality reads.
:::


Further details of the implementation of each flag can be found by reading the [Process and Parameters](#process-target) section below.  

All flags are tunable such that their parameters can be configured to a variety of use cases and sequencing methods.

(process-target)=
## Processes and Parameters

For an input VCF, `hairpin2` will iterate over the records therein and analyse each record via a series of interdependent scientific steps.

This section details each process and its relevant parameters (if any) and outcomes, in execution order. The connection/interdependence between steps is also described. The parameters are also shown in their relevant TOML table as written in a hairpin2 TOML config at the top of each section.

For each variant examined, `haripin2` determines the mutation type (SUB, INS, DEL), fetches all reads covering the mutant position from the alignments, and then walks through the steps described in this section. A reference table is provided below for a high level overview to refer back to while reading about each process in more detail.  

### Process Relationships Reference

:::{list-table} **Read Taggers**
:header-rows: 1

*   - Name
    - Tests Reads With
    - Excludes Reads With
    - Adds Tag
*   - mark-support
    - Any
    - None
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
*   - mark-overlap
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - None
    - {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`
*   - mark-low-qual
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - None
    - {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`
*   - mark-duplicates
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`
    - {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`
:::
:::{list-table} **Variant Flaggers**
:header-rows: 1

*   - Flag
    - Tests Reads With
    - Excludes Reads With
    - Examines Tag
*   - LQF
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`
    - {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`, {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`
*   - DVF
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`, {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`
    - {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`
*   - ALF
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`, {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`, {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`
    - None
*   - ADF
    - {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`
    - {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`, {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`, {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`
    - None
:::

### Read Taggers
---------

These processes execute prior to flaggers, examining the reads which cover the variant position. They apply extended tags to them such that downstream processes can filter/exclude/include reads from their analyses as needed. Note that tags applied to reads are internal, and are not written back to the alignment file.

#### mark-support

`mark-support` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagSupportingReads` to tag supporting reads for the variant in question.

Reads found to support the variant are marked with the tag {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG`

mark-support has no dependence on other steps

#### mark-overlap

`mark-overlap` takes no parameters. This process uses {py:class}`~hairpin2.sci_funcs.TagFragmentOverlapReads` to tag overlapping members of fragments (read pairs).

The second read found for a given qname is marked with the tag {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG`

mark-overlap only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support) - i.e. those that supporting the variant.

#### mark-low-qual
```toml
[params.mark-low-qual]
min_avg_clip_quality = 35
min_mapping_quality = 11
min_base_quality = 25
```

`mark-low-qual` uses {py:class}`~hairpin2.sci_funcs.TagLowQualReads` to mark reads that fail quality assessment.

Reads found to be of insufficient quality are marked with the tag {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG`. This step is essential to the LQF flag.

mark-low-qual only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. those that support the variant.

:::{list-table} **Parameters**
:header-rows: 1
:width: 100
:widths: 30 70
*   - Name
    - Description
*   - min_avg_clip_quality
    - Minimum mean quality of aligned bases for a read which contains soft clipping.
*   - min_mapping_quality
    - Minimum mapping quality (aka mapq) for a read.
*   - min_base_quality
    - Minimum base quality for any position covering the variant in a read.
:::

#### mark-duplicates
```toml
[params.mark-duplicates]
duplication_window_size = 6
```

(dup-explain-target)=
:::{admonition} Hidden Stutter Duplicates
In regions of low complexity, homopolymer tracts and secondary structure can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this variation, these duplicates both evade dupmarking and give rise to spurious variants when calling.
:::

`mark-duplicates` uses {py:class}`~hairpin2.sci_funcs.TagStutterDuplicateReads` to mark reads which appear to be hidden PCR duplicates.
The algorithm is simply based on assessing the variation of endpoints amongst the tested reads. If a read start/end and it's mate start/end
are within `duplication_window_size` of the closest neighbouring read's start/end and mate start/end then those reads are considered to
be duplicates of each other.

Reads found to be duplicates are marked with the tag {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`, excepting the highest quality read of the suspected duplicates. This step is essential to the DVF flag.

mark-duplicates only operates on reads marked with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. those that support the variant.

:::{list-table} **Parameters**
:header-rows: 1
:width: 100
:widths: 30 70
*   - Name
    - Description
*   - duplication_window_size
    - Alignment start/end window within which groups of fragments may be considered to be duplicates of one another.
:::


### Variant Flaggers
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


#### LQF
```toml
[params.LQF]
read_loss_threshold = 0.99
min_pass_reads = 2
nsamples_threshold = 0
```

`LQF` uses {py:class}`~hairpin2.sci_funcs.ProportionBasedTest` to flag variants which are largely supported by low quality reads.  

The test is performed by examining the proportion of reads supporting a variant marked with {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates) and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), compared to the total number of supporting reads (i.e. those tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support)).  

LQF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap) to avoid double counting.  

:::{list-table} **Parameters**
:header-rows: 1
:width: 100
:widths: 30 70

* - **Name**
  - **Description**
* - read_loss_threshold
  - Threshold of low-quality reads against total reads. If this threshold is exceeded, there are too many low-quality reads, and the variant is flagged.
* - min_pass_reads
  - Minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` and {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.
:::

:::{list-table} **INFO Conditions**
:header-rows: 1
:width: 100
:widths: 30 70

*   - Name
    - Description
*   - NO_READS
    - Check for suffcient reads available for testing, after excludes.
*   - THRESHOLD
    - Check for proportion of reads with stutter or low quality tags against the provided threshold.
*   - MIN_PASS
    - Check for sufficient reads without stutter or low quality tags against the provided minimum.
:::

In the case of a LQF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` may have `THRESHOLD`, `MIN_PASS`, or both. `NA` will only be seen with `NO_READS`. The INFO field for LQF additionally contains the proportion of reads tested found to have the {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`.  


#### DVF
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
:width: 100
:widths: 30 70
*  - Name
   - Description
*  - read_loss_threshold
   - threshold of duplicate reads against total reads.
    If this threshold is exceeded, there are too many duplicate reads,
    and the variant is flagged.
*  - min_pass_reads
   - minimum number of supporting reads without {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` necessary for a variant to pass the test. If there are fewer, the variant is flagged.
:::
:::{list-table} **INFO Conditions**
:header-rows: 1
:width: 100
:widths: 30 70

*   - Name
    - Description
*   - NO_READS
    - Check for suffcient reads available for testing, after excludes.
*   - THRESHOLD
    - Check for proportion of reads with the stutter tag against the provided threshold.
*   - MIN_PASS
    - Check for sufficient reads without the stutter tag against the provided minimum.
:::

In the case of a DVF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` may have `THRESHOLD`, `MIN_PASS`, or both. `NA` will only be seen with `NO_READS`. The INFO field for DVF additionally contains the proportion of reads tested found to have the {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG`.  

#### ALF
```toml
[params.ALF]
avg_AS_threshold = 0.93
```

`ALF` uses {py:class}`~hairpin2.sci_funcs.AlignmentScoreTest` to flag variants where the mean alignment score per base is less than `avg_AS_threshold`. The calcuation for each read is simply the value retrieved from the `AS` SAM tag divided by the read length. Users familiar with previous incarnations of LCMB variant post-processing should note that this is largely equivalent to `ASRD`.  

ALF operates only on reads tagged with {py:attr}`~hairpin2.const.Tags.SUPPORT_TAG` by [mark-support](#mark-support), i.e. only reads that support the given variant.

ALF excludes from testing any reads tagged with {py:attr}`~hairpin2.const.Tags.LOW_QUAL_TAG` by [mark-low-qual](#mark-low-qual), {py:attr}`~hairpin2.const.Tags.OVERLAP_TAG` by [mark-overlap](#mark-overlap), or {py:attr}`~hairpin2.const.Tags.STUTTER_DUP_TAG` by [mark-duplicates](#mark-duplicates).

:::{list-table} **Parameters**
:header-rows: 1
:width: 100
:widths: 30 70

*  - Name
   - Description
*  - avg_AS_threshold
   - threshold of mean alignment score per base. If this threshold is not reached, mean alignment score per base is too low, and the variant is flagged.
:::
:::{list-table} **INFO Conditions**
:header-rows: 1
:width: 100
:widths: 30 70

*   - Name
    - Description
*   - NO_READS
    - Check for suffcient reads available for testing, after excludes.
*   - INSUFFICIENT_AS_TAGS
    - Check for sufficient reads with AS tag to perform the test
*   - ON_THRESHOLD
    - Check for sufficient reads remaining without stutter or low quality tags
:::

In the case of a ALF `PASS`, all three condiitions will be noted in the INFO field, since the variant must necessarily have passed all three. A `FAIL` will only have `ON_THRESHOLD`. `NA` will be seen with one of `NO_READS` or `INSUFFICIENT_AS_TAGS`. The INFO field for ALF additionally contains the mean alignment score per base seen on reads tested for the variant.

#### ADF
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
:class: quote

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
:width: 100
:widths: 40 60
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

:::{list-table} **INFO Conditions**
:header-rows: 1
:width: 100
:widths: 30 70

*   - Name
    - Description
*   - NO_READS
    - Check for reads available for testing, after excludes.
*   - INSUFFICIENT_READS
    - Check for minimum reads available for testing on at least one strand, as defined by `low_n_supporting_reads_boundary`.
*   - EDGE_CLUSTERING
    - Check for fewer than `edge_clustering_threshold`% of reads expressing variant within `edge_definition` of alignment start, from both Option A and B in Ellis et al.
*   - ONE_STRAND_DISTRIB
    - Check for MAD and standard deviation condition from Option A in Ellis et. al., as defined by `min_MAD_one_strand` and `min_sd_one_strand`.
*   - BOTH_STRAND_DISTRIB_BOTH
    - Check for MAD and standard deviation condition from Option B as applied to both strands in Ellis et al., as defined by `min_MAD_both_strand_weak` and `min_sd_both_strand_weak`.
*   - BOTH_STRAND_DISTRIB_ONE
    - Check for MAD and standard devation condition from Option B as applied to one strand in Ellis et al., as defined by `min_MAD_both_strand_strong` and `min_sd_both_strand_strong`.
*   - MIN_NON_EDGE
    - Check for suffcient reads expressing variant away from read edge, as defined by `edge definition`.
:::

In the case of an ADF `PASS`/`FAIL`, all associated conditions passed/failed will be found in the INFO field. The specific conditions therein depend on which branch of the conditional was used. The INFO field contains the strand/s to help indicate the branch taken. `NA` will be seen with either NO_READS or INSUFFICIENT_READS only.  

## Writing a Config

hairpin2 configs can be written in either JSON or TOML format. We recommend [TOML](https://toml.io/en/), but both are fully supported. In either case, for standard usage a config needs a single top level key, `params`, which will then contain keys for each process by title, e.g. `mark-duplicates`, `ADF`, etc. Finally, under the process keys come the parameters. This format is pleasant to express in TOML:

```toml
[params.mark-duplicates]
duplication_window_size = 6


[params.DVF]
read_loss_threshold = 0.49
min_pass_reads = 2
nsamples_threshold = 0

# other processes... (n.b. toml allows comments)
```

if you prefer JSON:

```json
{
  "params": {
    "mark-low-qual": {
      "min_avg_clip_quality": 35,
      "min_mapping_quality": 11,
      "min_base_quality": 25
    },
    "LQF": {
      "read_loss_threshold": 0.99,
      "min_pass_reads": 2,
      "nsamples_threshold": 0
    }
  }
}
```

LCMB default parameters are provided in both JSON and TOML format in the `example-configs/` directory found in the hairpin2 [GitHub repository](https://github.com/cancerit/hairpin2).  

(understanding-decisions-target)=
## Understanding Decisions

hairpin2 stores detailed reasoning for each flagging decision in the INFO field of each tested variant. There will be VCF-spec standard key=value pairs for each flag tested. The key is the flag name, as would be applied to the FILTER column if the variant fails, e.g., `LQF`. For each flag, there will be as many key=value pairs as there are alts for the variant.  

The value side of the key=value pair is formatted as pipe (`|`) separated fields. The standard fields shared by all flags are as follows:

:::{list-table} **Value Fields**
:width: 100
:widths: 30 70
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
## Reproducibility & Parameter Distribution

hairpin2 uses only deterministic processes and so is reproducible within the standard computational bounds of that term. Singuarity and Docker definition files are provided for containerisation if needed.

hairpin2 stores the complete information necessary to distribute and reproduce a run directly in the VCF header. The following header keys are used

:::{list-table} **Header Keys**
:width: 100
:widths: 30 70
*   - hairpin2_version
    - stores version of hairpin2 used
*   - hairpin2_params
    - stores paramters and execution configuration
*   - hairpin2_samples
    - stores which samples in the VCF were selected for testing
:::

The parameters are stored in a reduced represenation and can be extracted from a VCF using [get-params](#get-params-qs-target) back into standard JSON format. It is not necessary to understand the way in which the parameters are encoded into the header (explained below) to use this functionality.  

The hairpin2_params key stores data in a JSON compressed with `zlib` and encoded to an ascii string using `base85`. The result is a string of ascii characters, much shorter than the original JSON, which does not invite or allow accidental editing of the parameters once written into the header (either by hand or by another tool). Since zlib compression includes an error-checking checksum, the string is guaranteed to transform back into the exact parameters encoded. Base85 encoding ensures that only VCF-safe characters are used  per the VCF format spec. 

## Assumptions & Limitations

hairpin2 is designed for paired-end data where alignment records have the `MC` tag and the complete CIGAR string is present in the `CIGAR` field (rather than the `CG:B,I` tag). If the `MC` tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. The tool can handle substitions, insertions, and deletions formatted per the VCF specification. At this time, the tool will not investigate mutations notated with angle brackets, e.g. `<DEL>`, complex mutations, or monomorphic reference. No further assumptions are made – other alignment tags and VCF fields are used, however they are mandatory per the relevant format specifications. If these requirements are limiting and you need the tool to be extended in some way, please request it.

## Further Development

We are very open to making updates and improvements to the tool to support a wide variety of use cases. If you have a request, please get in touch via the [GitHub](https://github.com/cancerit/hairpin2).

