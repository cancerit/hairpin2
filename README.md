# hairpin2
`hairpin2` – read-aware artefactual variant flagging

`hairpin2` is designed to flag variants that are likely artefactual via series of tests performed upon the read data associated with each variant. Initially, it was concieved to flag possible cruciform artefacts for LCM sequence data, but the concept has been extended and can detect a variety of potentially spurious variants (including indels). The tool operates on a VCF file containing one or more samples, and alignment files for all samples to be tested.

Given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with `ADF` if variants have anomalous distributions indicating that they are likely to be artefactual, `ALF` if relevant reads have lower median alignment score per base than a specified threshold, and `DVF` if variants appear to be the result of PCR error.

The `ADF ` filter is an implementation of the artifact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads.
The `ALF` filter indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` filter – artefacts with anomalous distributions often cause a marked decrease in alignment score.
The `DVF` filter is a naive but effective algorithm for detecting variants which are the result of PCR error - in regions of low complexity, short repeats and homopolymer tract can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling.

All filters are tunable such that their parameters can be configured to a variety of use cases and sequencing methods

### DEPENDENCIES
* `Python >= 3.12`

further dependencies (pysam, pydantic, and optionally pytest and pytest-cov) are detailed in `pyproject.toml`, and will be downloaded automatically if following the recommend install process

### INSTALLATION
The easiest end-user approach is to install into a virtual environment:
```
python -m venv .env
source .env/bin/activate
pip install .
hairpin -h
```

for development, substitute:
```
pip install -e ".[dev]"
```
to install testing dependencies


### ASSUMPTIONS & LIMITATIONS
`hairpin2` is designed for paired data where alignment records have the `MC` tag and the complete CIGAR string is present in the `CIGAR` field (rather than the `CG:B,I` tag). If the `MC` tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. The tool can handle substitions, insertions, and deletions formatted per the VCF specification. At this time, the tool will not investigate mutations notated with angle brackets, e.g. `<DEL>`, complex mutations, or monomorphic reference. No further assumptions are made – other alignment tags and VCF fields are used, however they are mandatory per the relevant format specifications. If these requirements are limiting and you need the tool to be extended in some way, please request it.

### USAGE
```
Usage: hairpin2 [-h, --help] [OPTIONS] VCF_IN ALIGNMENTS...

  read-aware artefactual variant flagging algorithms. Flag variants in VCF
  using statistics calculated from supporting reads found in ALIGNMENTS, and
  emit the flagged VCF to stdout.

Options:
  -v, --version                   Show the version and exit.
  -h, --help                      show helptext
  -hh, --help-all                 Show further help including config override
                                  options
  -r, --cram-reference FILEPATH   path to FASTA format CRAM reference,
                                  overrides $REF_PATH and UR tags - ignored if
                                  --format is not CRAM
  -m, --name-mapping S:SM S:SM...
                                  key to map samples in a multisample VCF to
                                  alignment/s provided to -a. Uses VCF sample
                                  names from VCF header and alignment SM tags.
                                  With multiple alignments to -a, accepts a
                                  space separated list of sample:SM pairs.
                                  When only a single alignment provided, also
                                  accepts a comma separated string of one or
                                  more possible sample-of-interest names like
                                  TUMOR,TUMOUR
  -c, --config FILEPATH           path to config JSON from which filter
                                  paramters will be loaded - can be overridden
                                  by extended arguments provided at runtime
  -o, --output_config FILEPATH    log filter paramaters from run as a config
                                  JSON file
  --progess / --no-progess        display progress bar on stderr during run

read validation config overrides:
    --min-clip-quality INT        discard reads with mean base quality of
                                  aligned bases below this value, if they have
                                  soft-clipped bases  [default: (35);
                                  0<=x<=60]
    --min-mapping-quality INT     discard reads with mapping quality below
                                  this value  [default: (11); 0<=x<=93]
    --min-base-quality INT        discard reads with base quality below this
                                  value at variant position  [default: (25);
                                  x>=-1]
  DVF config overrides:
    --duplication-window-size INT
                                  inclusive maximum window size, in number of
                                  bases to use when detecting PCR duplicates.
                                  -1 will disable duplicate detection
                                  [default: (6); x>=0]

ALF config overrides:
    --al-filter-threshold FLOAT   ALF; threshold for median of read alignment
                                  score per base of all relevant reads, at and
                                  below which a variant is flagged as ALF
                                  [default: (0.93); 0<=x<=93]

ADF config overrides:
    --edge-definition FLOAT       ADF; percentage of a read that is considered
                                  to be "the edge" for the purposes of
                                  assessing variant location distribution
                                  [default: (0.15); 0.0<=x<=0.99]
    --edge-fraction FLOAT         ADF; percentage of variants must occur
                                  within EDGE_FRACTION of read edges to mark
                                  ADF flag  [default: (0.9); 0.0<=x<=0.99]
    --min-mad-one-strand INT      ADF; min range of distances between variant
                                  position and read start for valid reads when
                                  only one strand has sufficient valid reads
                                  for testing  [default: (0); x>=0]
    --min-sd-one-strand FLOAT     ADF; min stdev of variant position and read
                                  start for valid reads when only one strand
                                  has sufficient valid reads for testing
                                  [default: (4.0); x>=0.0]
    --min-mad-both-strand-weak INT
                                  ADF; min range of distances between variant
                                  position and read start for valid reads when
                                  both strands have sufficient valid reads for
                                  testing AND -sbsw is true  [default: (2);
                                  x>=0]
    --min-sd-both-strand-weak FLOAT
                                  ADF; min stdev of variant position and read
                                  start for valid reads when both strands have
                                  sufficient valid reads for testing AND -mbsw
                                  is true- default: 2, range: 0-, exclusive
                                  [default: (2.0); x>=0.0]
    --min-mad-both-strand-strong INT
                                  ADF; min range of distances between variant
                                  position and read start for valid reads when
                                  both strands have sufficient valid reads for
                                  testing AND -sbss is true  [default:
                                  (ParamConstraint(default=1,
                                  range=MinMax(min=0, max=None))); x>=0]
    --min-sd-both-strand-strong FLOAT
                                  ADF; min stdev of variant position and read
                                  start for valid reads when both strands have
                                  sufficient valid reads for testing AND -mbss
                                  is true  [default: (10.0); x>=0.0]
    --min-reads INT               ADF; number of reads at and below which the
                                  hairpin filtering logic considers a strand
                                  to have insufficient reads for testing
                                  [x>=1]
```

Parameters are hopefully mostly clear from the helptext, but some warrant further explanation:

- `--name-mapping` – When using multisample VCFS, hairpin2 compares VCF sample names found in the VCF header to SM tags in alignments to match samples of interest to the correct alignment. If these IDs are different between the VCF and alignments, you'll need to provide a key. If there are multiple samples of interest in the VCF, and therefore multiple alignments, you will need to provide a key for each pair - e.g. `-m sample1:SM1 sample2:SM2 ...`. If there is only one alignment, then you need only indicate which VCF sample is the sample of interest, e.g. `-m TUMOR`. As a convenience for high throughput workflows, when there is only one alignment you may also provide a comma separated string of possible names for the sample of interest, e.g. `-m TUMOR,TUMOUR`. Assuming there is one and only one match in the VCF, the tool will match the alignment to that sample.
- `--al-filter-threshold` – the default value of 0.93 was arrived at by trial and error – since different aligners/platforms calculate alignment score differently, you may want to modify this value appropriately. In the predecessor to this tool, `additionalBAMStatistics`, this value was known as `ASRD` and the default was set at 0.87.
- `--duplication-window-size` – long homopolymer tracts can cause stuttering, where a PCR duplicate will have, for example, an additional A in a tract of As. These reads will align a base or two earlier on the reference genome than they should. As a result pcr duplicate flag machinery fails and they are not flagged as duplicates. `hairpin2` will attempt to filter out these duplicates, and duplication window size is then the maximum +- position to use during duplicate detection.

The parameters available for the ADF flag are probably best understood by reading the 'test()' method implementation of that filter.

The tool tests records in a VCF file and applies filter flags as appropriate. Reasoning for decisions is recorded in the INFO field of the VCF records, in the form `<FILTER>=<alt>|<True/False>|<code>|...`. The codes indicate the reason on which the decision was made, and are as follows:

> **0** – passed/failed on condition 60A(i) of Ellis et al. (`ADF` only)  
> **1** – passed/failed on condition 60B(i) of Ellis et al. (`ADF` only)  
> **2** – passed/failed on filter threshold (`ALF` only)  
> **3** – insufficient appropriate reads to support calling flag – this covers a lot of possiblities, if more granularity is desired, please request it  
> **4** – no samples have non 0,0 genotype for the record  

The basic procedure of this implementation is as follows:  
>   For each record in the VCF, test every alt for that record as follows:  
>   1. for samples exhibiting the mutation, retrieve reads covering the region
>   2. test each read for validity for use in testing (i.e. base quality, do they express the correct alt, and so on)
>   3. performing statistical analysis on read context
>   4. pass or fail the record for each filter and log to `INFO`


### TESTING
A test suite has been provided with the algorithm implementation. To run these tests run `pytest .` from within the install directory. `hairpin2` must have been installed from that same directory, and be available on path (for example in a virtual environment). The tests can be found in the `test` directory. The focus is on testing all nodes (statements) and edges (paths between statements) for key scientific logic, rather than trying every possible input combination, in addition to testing critical boundary conditions. Since all input possibilities will pass via the same network/graph, if we prove that each part of that network functions correctly, we can be more confident that the program functions as it claims to. The tests are simple, and, once you have read them, it should be very easy to add your own further tests should you feel the need to confirm any behaviour.


### TODO
- automated regression testing
- stricter config validation, most likely with pydantic, to catch misformatted configs earlier
- improve documentation - describe filters in individual sections, beyond replicating helptext


### LICENCE
```
hairpin2

Copyright (C) 2024, 2025 Genome Research Ltd.

Author: Alex Byrne <ab63@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
