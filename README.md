# hairpin2

`hairpin2` – CLI implementation of the hairpin detection algorithm concieved by [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). 

For paired data, given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with **HPF** if they are suspected cruciform artefacts, and **ALF** if relevant reads have lower median alignment score per base than a specified threshold. The **ALF** filter indicates poor signal-to-noise, and provides additional confidence in the **HPF** filter – cruciform artefacts usually cause a marked decrease in alignment score. The **ALF** flag also may appear on variants without **HPF**, often indicating other artefacts associated with poor signal-to-noise.


### DEPENDENCIES

* Python >= 3.10 – required
* pysam >= 0.22.1 – installed automatically during install process (tested with 0.22.1 only)

### INSTALLATION


Within a virtual environment:
```
python -m venv .env
source .env/bin/activate
pip install .
hairpin -h
```

For system-wide access:
```
export INST_PATH=/path/to/install/location/
mkdir -p $INST_PATH
pip install . --target $INST_PATH
export PATH=${PATH}:${INST_PATH}/bin
hairpin -h
```

### ASSUMPTIONS

`hairpin2` is designed for paired data where BAM records have the **MC** tag. If this tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. The tool expects data specifically in the VCF and BAM formats; support for a wider variety of formats could be implemented if desired. No further assumptions are made – other BAM tags and VCF fields are used, however they are mandatory per the format specification.


### USAGE

```
usage: hairpin2 [-h] [-v] [-i VCF_IN] [-o VCF_OUT] [-b BAMS [BAMS ...]] [-al AL_FILTER_THRESHOLD] [-mc MIN_CLIP_QUALITY] [-mq MIN_MAPPING_QUALITY] [-mb MIN_BASE_QUALITY]
                [-ms MAX_READ_SPAN] [-pf POSITION_FRACTION] [-m VCF:BAM [VCF:BAM ...]] [-ji INPUT_JSON] [-jo OUTPUT_JSON]

cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6)

info:
  -h, --help            show this help message and exit
  -v, --version         print version

basic:
  -i VCF_IN, --vcf-in VCF_IN
                        path to input VCF
  -o VCF_OUT, --vcf-out VCF_OUT
                        path to write output VCF
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        list of paths to BAMs for samples in input VCF, whitespace separated

extended:
  -al AL_FILTER_THRESHOLD, --al-filter-threshold AL_FILTER_THRESHOLD
                        threshhold for median of read alignment score per base of all relevant reads, below which a variant is flagged as ALF - default: 0.93
  -mc MIN_CLIP_QUALITY, --min-clip-quality MIN_CLIP_QUALITY
                        discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases - default: 35
  -mq MIN_MAPPING_QUALITY, --min-mapping-quality MIN_MAPPING_QUALITY
                        discard reads with mapping quality below this value - default: 11
  -mb MIN_BASE_QUALITY, --min-base-quality MIN_BASE_QUALITY
                        discard reads with base quality at variant position below this value - default: 25
  -ms MAX_READ_SPAN, --max-read-span MAX_READ_SPAN
                        maximum +- position to use when detecting PCR duplicates - default: 6
  -pf POSITION_FRACTION, --position-fraction POSITION_FRACTION
                        >90% of variant must occur within POSITION_FRACTION of read edges to allow HPF flag - default: 0.15

procedural:
  -m VCF:BAM [VCF:BAM ...], --name-mapping VCF:BAM [VCF:BAM ...]
                        map VCF sample names to BAM SM tags; useful if they differ
  -ji INPUT_JSON, --input-json INPUT_JSON
                        path to JSON of input parameters; overridden by arguments provided on command line
  -jo OUTPUT_JSON, --output-json OUTPUT_JSON
                        log input arguments to JSON
```

Parameters are hopefully mostly clear from the helptext, but some warrant further explanation:

> `--al-filter-threshold` – the default value of 0.93 was arrived at by trial and error – since different aligners/platforms calculate alignment score differently, you may want to modify this value appropriately.  
> `--max-read-span`  – long homopolymer tracts can cause stuttering, where a PCR duplicate will have, for example, an additional A in a tract of As. These reads will align a base or two earlier on the reference genome than they should. As a result pcr duplicate flag machinery fails and they are not flagged as duplicates. `MAX_READ_SPAN` is then the maximum +- position to use when detecting PCR duplicates.  
> `--position-fraction` – cruciform artefacts usually contain segments that do not align to the reference genome, resulting in the segment being soft-clipped. The subsequent aligned portion will then contain false variants, which arise from the artefact. These false variants appear with anomalous regularity at alignment boundaries – unlike true variants. If, for a given variant, more than 90% of the variant bases are within `POSITION_FRACTION` of read edges, allow for calling **HPF** flag.


### DETAILS

The tool tests records in a VCF file and applies the **HPF** and **ALF** filter flags as appropriate. Reasoning for decisions is recorded in the INFO field of the VCF records, in the form `HPF=<alt>|<code>` and `ALF=<alt>|<code>|<median AS score>`. The codes are as follows:  

> **0** – passed/failed on condition 60A(i) of Ellis et al. (HPF only)  
> **1** – passed/failed on condition 60B(i) of Ellis et al. (HPF only)  
> **2** – passed/failed on filter threshold (ALF only)  
> **3** – insufficient appropriate reads to support calling flag (pass only)   (This covers a lot of possiblities, if more granularity is desired, please request it)  
> **4** – no samples have non 0,0 genotype for the record (pass only)
  

The basic procedure of this implementation is as follows:  
>   For each record in the VCF, test every alt for that record by:  
>   1. retrieving reads from samples exhibiting the mutations
>   2. testing each read for validity for use in hairpin testing (i.e. base quality, do they express the correct alt, and so on)
>   3. performing statistical analysis on aggregates of the position of the mutation relative to the start and end of the aligned portion of the reads
>   4. on the results of the statistical analysis, pass or fail the record for the filters **ALF** and **HPF**, and log a code and relevant info to the **INFO** field indicating the reason for the decision