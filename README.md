# hairpin2

Maintainable, transparent, implementation of the hairpin detection and flagging algorithm concieved by Mathijs' Sanders. Implemented by Peter Campbell and Alex Byrne

### REQUIREMENTS

* Python >= 3.10

### INSTALLATION


Within a virtual environment:
```
python -m venv .env
source .env/bin/activate
pip install .
hairpin2 -h
```

For system-wide access:
```
export INST_PATH=/path/to/install/location/
mkdir -p $INST_PATH
pip install . --target $INST_PATH
export PATH=${PATH}:${INST_PATH}/bin
hairpin2 -h
```

### DETAILS

```
usage: hairpin2 [-h] [-v] -i VCF_IN -o VCF_OUT -b BAMS [BAMS ...] [-cq CLIP_QUALITY_CUTOFF] [-mq MIN_MAPPING_QUALITY] 
  [-mb MIN_BASE_QUALITY] [-ms MAX_READ_SPAN] [-al AL_FILTER_THRESHOLD] [-c9 CENT90_THRESHOLD] [-j JSON_PATH]

  info:
    -h, --help            show this help message and exit
    -v, --version         print version

  required:
    -i VCF_IN, --vcf-in VCF_IN
                path to input vcf
    -o VCF_OUT, --vcf-out VCF_OUT
                path to vcf out
    -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                list of paths to bams for samples in input vcf, whitespace separated

  options:
    -cq CLIP_QUALITY_CUTOFF, --clip-quality-cutoff CLIP_QUALITY_CUTOFF
                default: 35
    -mq MIN_MAPPING_QUALITY, --min-mapping-quality MIN_MAPPING_QUALITY
                default: 11
    -mb MIN_BASE_QUALITY, --min-base-quality MIN_BASE_QUALITY
                default: 25
    -ms MAX_READ_SPAN, --max-read-span MAX_READ_SPAN
                default: 6
    -al AL_FILTER_THRESHOLD, --al-filter-threshold AL_FILTER_THRESHOLD
                default: 0.93
    -c9 CENT90_THRESHOLD, --cent90-threshold CENT90_THRESHOLD
                default: 0.15
    -j JSON_PATH, --json-log JSON_PATH
                log input parameters/arguments to JSON
```

The tool tests records in a VCF file and applies the HPF, indicating a hairpin, and ALF, flags as appropriate. It records reasoning for its decisions in the INFO field of the VCF records, in the form HPF=<alt>|<code> and ALF=<alt>|<code>|<average AS score>

The codes are as follows

0 - passed/failed on condition 60A(i) of Ellis et al. (HPF only)

1 - passed/failed on condition 60B(i) of Ellis et al. (HPF only)

2 - passed/failed on filter threshold (ALF only)

3 - insufficient appropriate reads to support calling flag (This covers a lot of possiblities, if more granularity is desired, please request it)

4 - no samples have non 0,0 genotype for the record

The basic procedure of this implementation is as follows:
> For each record in the VCF, test every alt for that record by:
> * retrieving reads from samples exhibiting the mutations, then
> * testing each read for validity for use in hairpin testing (i.e. base quality, do they express the correct alt, and so on), then
> * performing statistical analysis on aggregates of the position of the mutatation relative to the start and end of the aligned portion of the reads, then
> * on the results of the statistical analysis, pass or fail the record for the filters ALF and HPF, and log a code and relevant info to the INFO field indicating the reason for the decision
