# hairpin2

`hairpin2` – CLI implementation of the hairpin detection algorithm concieved by [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6).

`hairpin2` is designed to flag variants as possible cruciform artefacts. It operates on a VCF file containing one or more samples, and alignment files for all samples to be tested.

Given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with `HPF` if they are suspected cruciform artefacts, and `ALF` if relevant reads have lower median alignment score per base than a specified threshold.

The `ALF` filter indicates poor signal-to-noise, and provides additional confidence in the `HPF` filter – cruciform artefacts usually cause a marked decrease in alignment score. The `ALF` flag also may appear on variants without `HPF`, often indicating other artefacts associated with poor signal-to-noise.


### DEPENDENCIES

* Python >= 3.10 – required
* pysam >= 0.22.1 – installed automatically during install process (tested with 0.22.1 only)
* pytest >= 0.8.2.2 - optional, only necessary to run tests

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

`hairpin2` is designed for paired data where alignment records have the `MC` tag and the complete CIGAR string is present in the `CIGAR` field (rather than the `CG:B,I` tag). If the `MC` tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. No further assumptions are made – other alignment tags and VCF fields are used, however they are mandatory per the relevant format specifications.


### USAGE

```
usage: hairpin2 [-h] [-v] -i VCF_IN -o VCF_OUT -a ALIGNMENTS [ALIGNMENTS ...]
                -f {s,b,c} [-al AL_FILTER_THRESHOLD] [-mc MIN_CLIP_QUALITY]
                [-mq MIN_MAPPING_QUALITY] [-mb MIN_BASE_QUALITY]
                [-ms MAX_READ_SPAN] [-pf POSITION_FRACTION]
                [-r CRAM_REFERENCE] [-m VCF:aln [VCF:aln ...]]
                [-ji INPUT_JSON] [-jo OUTPUT_JSON]

cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI:
10.1038/s41596-020-00437-6). See README for further explanation of parameters.

info:
  -h, --help            show this help message and exit
  -v, --version         print version

mandatory:
  -i VCF_IN, --vcf-in VCF_IN
                        path to input VCF
  -o VCF_OUT, --vcf-out VCF_OUT
                        path to write output VCF
  -a ALIGNMENTS [ALIGNMENTS ...], --alignments ALIGNMENTS [ALIGNMENTS ...]
                        list of paths to (S/B/CR)AMs (indicated by --format)
                        for samples in input VCF, whitespace separated -
                        (s/b/cr)ai expected in same directories
  -f {s,b,c}, --format {s,b,c}
                        format of alignment files; s indicates SAM, b
                        indicates BAM, and c indicates CRAM

read validation:
  -mc MIN_CLIP_QUALITY, --min-clip-quality MIN_CLIP_QUALITY
                        discard reads with mean base quality of aligned bases
                        below this value, if they have soft-clipped bases -
                        default: 35, range: 0-93, exclusive
  -mq MIN_MAPPING_QUALITY, --min-mapping-quality MIN_MAPPING_QUALITY
                        discard reads with mapping quality below this value -
                        default: 11, range: 0-60, exclusive
  -mb MIN_BASE_QUALITY, --min-base-quality MIN_BASE_QUALITY
                        discard reads with base quality at variant position
                        below this value - default: 25, range: 0-93, exclusive
  -ms MAX_READ_SPAN, --max-read-span MAX_READ_SPAN
                        maximum +- position to use when detecting PCR
                        duplicates. -1 will disable duplicate detection -
                        default: 6, range: -1-, inclusive

filter conditions:
  -al AL_FILTER_THRESHOLD, --al-filter-threshold AL_FILTER_THRESHOLD
                        ALF; threshold for median of read alignment score per
                        base of all relevant reads, at and below which a
                        variant is flagged as ALF - default: 0.93, range: 0-,
                        inclusive
  -ed EDGE_DEFINITION, --edge-definition EDGE_DEFINITION
                        HPF; percentage of a read that is considered to be
                        "the edge" for the purposes of assessing variant
                        location distribution - default: 0.15, range: 0-0.99,
                        inclusive
  -ef EDGE_FRACTION, --edge-fraction EDGE_FRACTION
                        HPF; percentage of variants must occur within
                        EDGE_FRACTION of read edges to allow HPF flag -
                        default: 0.15, range: 0-0.99, exclusive
  -mos MIN_MAD_ONE_STRAND, --min-MAD-one-strand MIN_MAD_ONE_STRAND
                        HPF; min range of distances between variant position
                        and read start for valid reads when only one strand
                        has sufficient valid reads for testing - default: 0,
                        range: 0-, exclusive
  -sos MIN_SD_ONE_STRAND, --min-sd-one-strand MIN_SD_ONE_STRAND
                        HPF; min stdev of variant position and read start for
                        valid reads when only one strand has sufficient valid
                        reads for testing - default: 4, range: 0-, exclusive
  -mbsw MIN_MAD_BOTH_STRAND_WEAK, --min-MAD-both-strand-weak MIN_MAD_BOTH_STRAND_WEAK
                        HPF; min range of distances between variant position
                        and read start for valid reads when both strands have
                        sufficient valid reads for testing AND -sbsw is true -
                        default: 2, range: 0-, exclusive
  -sbsw MIN_SD_BOTH_STRAND_WEAK, --min-sd-both-strand-weak MIN_SD_BOTH_STRAND_WEAK
                        HPF; min stdev of variant position and read start for
                        valid reads when both strands have sufficient valid
                        reads for testing AND -mbsw is true- default: 2,
                        range: 0-, exclusive
  -mbss MIN_MAD_BOTH_STRAND_STRONG, --min-mad-both-strand-strong MIN_MAD_BOTH_STRAND_STRONG
                        HPF; min range of distances between variant position
                        and read start for valid reads when both strands have
                        sufficient valid reads for testing AND -sbss is true -
                        default: 1, range: 0-, exclusive
  -sbss MIN_SD_BOTH_STRAND_STRONG, --min-sd-both-strand-strong MIN_SD_BOTH_STRAND_STRONG
                        HPF; min stdev of variant position and read start for
                        valid reads when both strands have sufficient valid
                        reads for testing AND -mbss is true - default: 10,
                        range: 0-, exclusive
  -mr MIN_READS, --min-reads MIN_READS
                        HPF; number of reads at and below which the hairpin
                        filtering logic considers a strand to have
                        insufficient reads for testing - default: 1, range: 0-

procedural:
  -r CRAM_REFERENCE, --cram-reference CRAM_REFERENCE
                        path to FASTA format CRAM reference, overrides
                        $REF_PATH and UR tags - ignored if --format is not
                        CRAM
  -m VCF:aln [VCF:aln ...], --name-mapping VCF:aln [VCF:aln ...]
                        map VCF sample names to alignment SM tags; useful if
                        they differ
  -ji INPUT_JSON, --input-json INPUT_JSON
                        path to JSON of input parameters, from which extended
                        arguments will be loaded - overridden by arguments
                        provided on command line
  -jo OUTPUT_JSON, --output-json OUTPUT_JSON
                        log input arguments to JSON
```

Parameters are hopefully mostly clear from the helptext, but some warrant further explanation:

- --name-mapping – some variant callers, for example caveman, output sample names such as "TUMOUR" in VCF header columns. hairpin2 uses these column names to match to BAM samples via the SM tag - if these fields do not match, you'll need to provide a mapping here, for example "TUMOR:PD3738..."
- --al-filter-threshold – the default value of 0.93 was arrived at by trial and error – since different aligners/platforms calculate alignment score differently, you may want to modify this value appropriately. In "Mathijs' Scripts", the default was set at 0.87 for filtering on ASRD.
- --max-read-span – long homopolymer tracts can cause stuttering, where a PCR duplicate will have, for example, an additional A in a tract of As. These reads will align a base or two earlier on the reference genome than they should. As a result pcr duplicate flag machinery fails and they are not flagged as duplicates. `hairpin2` will attempt to filter out these duplicates, and MAX_READ_SPAN is then the maximum +- position to use during duplicate detection.
- --edge-fraction – cruciform artefacts usually contain segments that do not align to the reference genome, resulting in the segment being soft-clipped. The subsequent aligned portion will then contain false variants, which arise from the artefact. These false variants appear with anomalous regularity at alignment boundaries – unlike true variants. If, for a given variant, more than 90% of the variant bases are within EDGE_FRACTION of read edges, allow for calling HPF flag.

The parameters available for the HPF flag are probably best understood by reading the implementation of the function `is_variant_HP()` in `hairpin2/main.py`.

The tool tests records in a VCF file and applies the `HPF` and `ALF` filter flags as appropriate. Reasoning for decisions is recorded in the INFO field of the VCF records, in the form `HPF=<alt>|<True/False>|<code>` and `ALF=<alt>|<True/False>|<code>|<median AS score>`. The codes are as follows:  

> **0** – passed/failed on condition 60A(i) of Ellis et al. (`HPF` only)  
> **1** – passed/failed on condition 60B(i) of Ellis et al. (`HPF` only)  
> **2** – passed/failed on filter threshold (`ALF` only)  
> **3** – insufficient appropriate reads to support calling flag – this covers a lot of possiblities, if more granularity is desired, please request it  
> **4** – no samples have non 0,0 genotype for the record  

The basic procedure of this implementation is as follows:  
>   For each record in the VCF, test every alt for that record as follows:  
>   1. for samples exhibiting the mutation, retrieve reads covering the region
>   2. test each read for validity for use in hairpin testing (i.e. base quality, do they express the correct alt, and so on)
>   3. performing statistical analysis on aggregates of the position of the mutation relative to the start and end of the aligned portion of the reads
>   4. on the results of the statistical analysis, pass or fail the record for the filters `ALF` and `HPF`, and log a code and relevant info to the `INFO` field indicating the reason for the decision  

The code has been written with the intention of clarity and extensibility – again, further understanding may be achieved by reading `hairpin2/main.py`.


### TESTING

A test suite has been provided to prove the validity of the algorithm. To run these tests run `pytest -m "validate"` from within the install directory. `hairpin2` must have been installed from that same directory, and be available on path. The tests can be found in the `test` directory. They are simple, and, once you have read them, it should be very easy to add your own further tests should you want to confirm any behaviour.


### LICENCE

```
hairpin2

Copyright (C) 2024 Genome Research Ltd.

Author: Alex Byrne <ab63@sanger.ac.uk>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
