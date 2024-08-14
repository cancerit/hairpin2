### INTRODUCTION

`hairpin2` - CLI implementation of the hairpin detection and flagging algorithm concieved by [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). Implemented by Peter Campbell and Alex Byrne (primary contact for this tool - ab63). Code not yet public, but availabe on internal gitlab at https://gitlab.internal.sanger.ac.uk/casm/team78/hairpin-core

For paired data, given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with **HPF** if they are suspected cruciform artefacts, and **ALF** if they exhibit a lower median alignment score than a specified threshold. The **ALF** filter indicates poor signal-to-noise, and provides additional confidence in the **HPF** filter - cruciform artefacts usually cause a marked decrease in alignment score. The **ALF** flag also may appear on variants without **HPF**, often indicating other artefacts associated with poor signal-to-noise.

`hairpin2` should replace, as far as is possible, the tools known as "Mathijs' Scripts", "AdditionalBamStatistics", "Tim Butler's Scripts" and, unfortunately, probably many other names. It also supersedes `hairpin`, a stopgap version /of Mathijs' Scripts that relied on some of Mathijs' original code, and therefore was unreliable and error prone (though less so than the raw scripts themselves).
However, this incarnation is not a total replacement for Mathijs' Scripts at this time (and has changed in functionality since the stopgap tool, the original hairpin):

> """  
> Mathjis LCM filters includes the following steps:
> 1. Preselect: Filters the CaVEMan calls for “PASS” && “CLPM=0” && “ASMD>=140”
> 2. Hairpin Filtering
> 3. Filtering based on fragment numbers.  
>
> Which are split across the following steps: (As per his scripts)  
> - preselect
> - imitateANNOVAR
> - annotateBAMStatistics
> - additionalBAMStatistics
> - filtering  
>
> The `hairpin2` module replaces the “additionalBAMStatistics” and most of the “filtering” code. So [one may still need] to run the preselect and fragment based filter.  
> """  

Since the versions available of "Mathijs' Scripts" are many and varied, we cannot account for all differences/changes, but in general:
> - No more ambiguous/cryptic/unfixable errors - the tool should work on all appropriate data, and if it is unable to produce the expected output it will clearly inform the user (but see N.B. at end of this section)
> - Transparency - reasoning for flagging decisions logged in VCF
> - Single tool centrally maintained and versioned - for reproducibility/citing/distribution
> - Significant speedup (on testing data at least) – 50s runtime on 542-variant caveman VCF
> - The module adds **filter flags**, **HPF** and **ALF**, to a VCF. It **does not** output into separate files containing passed and failed positions
> - The module **does not** prefilter, or perform fragment filtering
> With regard to prefiltering - this is not performed by this module, as the filtering is not relevant to hairpin detection and should be performed separately. Filtering can be performed using the `vcfilter` or `bcftools` modules.  

**N.B.** this program is currently in an alpha/testing phase - it is available on the farm, but is likely to change, or have new features added, rapidly, per user responses. **It also may be broken in some way; if so please get in touch**. It is not currently publicly available - it will be made public as soon as it is out of this alpha phase.


### ACCESS

For local or VM use, see GitLab for install instructions.
For farm22 use, available as a module.
```
module avail hairpin2
module load <version>
```
 **N.B. do not confuse with the module `hairpin` - this is `hairpin2`**


### ASSUMPTIONS

`hairpin2` is designed for paired data where reads have the **MC** tag. If this tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. The tool expects data specifically in the VCF and BAM formats; support for a wider variety of formats could be implemented if desired.


### USAGE

```
usage: hairpin [-h] [-v] [-i VCF_IN] [-o VCF_OUT] [-b BAMS [BAMS ...]] [-al AL_FILTER_THRESHOLD] [-mc MIN_CLIP_QUALITY] [-mq MIN_MAPPING_QUALITY] [-mb MIN_BASE_QUALITY] [-ms MAX_READ_SPAN] [-pf POSITION_FRACTION]
               [-m VCF:BAM [VCF:BAM ...]] [-ji INPUT_JSON] [-jo OUTPUT_JSON]

info:
  -h, --help            show this help message and exit
  -v, --version         print version

basic:
  -i VCF_IN, --vcf-in VCF_IN
                        path to input vcf
  -o VCF_OUT, --vcf-out VCF_OUT
                        path to vcf out
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        list of paths to name-sorted bams for samples in input vcf, whitespace separated

extended:
  -al AL_FILTER_THRESHOLD, --al-filter-threshold AL_FILTER_THRESHOLD
                        threshhold for median of read alignment score per base of all relevant reads, below which a variant is flagged as **ALF** - default: 0.93
  -mc MIN_CLIP_QUALITY, --min-clip-quality MIN_CLIP_QUALITY
                        discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35
  -mq MIN_MAPPING_QUALITY, --min-mapping-quality MIN_MAPPING_QUALITY
                        discard reads with mapping quality below this value - default: 11
  -mb MIN_BASE_QUALITY, --min-base-quality MIN_BASE_QUALITY
                        discard reads with base quality at variant position below this value - default: 25
  -ms MAX_READ_SPAN, --max-read-span MAX_READ_SPAN
                        maximum +- position to use when detecting PCR duplicates - default: 6
  -pf POSITION_FRACTION, --position-fraction POSITION_FRACTION
                        >90% of variant reads variant must occur within [fraction] of start/end to allow **HPF** flag - default: 0.15

procedural:
  -m VCF:BAM [VCF:BAM ...], --name-mapping VCF:BAM [VCF:BAM ...]
                        map VCF sample names to BAM sample names; useful if they differ
  -ji INPUT_JSON, --input-json INPUT_JSON
                        path to JSON of input parameters; overridden by arguments provided on command line
  -jo OUTPUT_JSON, --output-json OUTPUT_JSON
```

**N.B.** the above usage block indicates the call for the tool is `hairpin2` - this is correct for local/vm installs, but for farm usage, for the time being, it is `hairpin2-alpha`

Parameters are hopefully mostly clear from the helptext, but some warrant further explanation:

> `--al-filter-threshold` - the default value of 0.93 was arrived at by trial and error - since different aligners/platforms calculate alignment score differently, you may want to modify this value appropriately.  
> `--max-read-span`  - Long homopolymer tracts can cause stuttering, where a PCR duplicate will have, for example, an additional A in a tract of As. These reads will align a base or two earlier on the reference genome than they should. As a result pcr duplicate flag machinery fails and they are not flagged as duplicates. max-read-span is then the maximum +- position to use when detecting PCR duplicates.  
> `--position-fraction` - cruciform artefacts usually contain segments that align beause the segment is not in ref genome, and so the segment is soft clipped – this pushes the false variants associated with the arterfact to edges of the reads; unlike true variants. If more than 90% of the reads are within that first/last fraction, allow for calling **HPF** flag



### DETAILS

The tool tests records in a VCF file and applies the **HPF**, indicating a hairpin/cruciform artefact, and **ALF** filter flags as appropriate. It records reasoning for its decisions in the INFO field of the VCF records, in the form `HPF=<alt>|<code>` and `ALF=<alt>|<code>|<median AS score>`. The codes are as follows:  

> **0** - passed/failed on condition 60A(i) of Ellis et al. (HPF only)  
> **1** - passed/failed on condition 60B(i) of Ellis et al. (HPF only)  
> **2** - passed/failed on filter threshold (ALF only)  
> **3** - insufficient appropriate reads to support calling flag (pass only)   (This covers a lot of possiblities, if more granularity is desired, please request it)  
> **4** - no samples have non 0,0 genotype for the record (pass only)
  

The basic procedure of this implementation is as follows:  
>   For each record in the VCF, test every alt for that record by:  
>   1. retrieving reads from samples exhibiting the mutations
>   2. testing each read for validity for use in hairpin testing (i.e. base quality, do they express the correct alt, and so on)
>   3. performing statistical analysis on aggregates of the position of the mutatation relative to the start and end of the aligned portion of the reads
>   4. on the results of the statistical analysis, pass or fail the record for the filters **ALF** and **HPF**, and log a code and relevant info to the **INFO** field indicating the reason for the decision