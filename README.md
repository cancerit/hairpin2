# hairpin2
`hairpin2` – read-aware artefactual variant flagging

`hairpin2` is designed to flag variants that are likely artefactual via a series of tests performed upon the read data associated with each variant. Initially, it was concieved to flag possible cruciform artefacts for LCM sequence data, but the concept has been extended and can detect a variety of potentially spurious variants (including indels). The tool operates on a VCF file containing one or more samples, and alignment files for all samples to be tested.

Given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with `ADF` if variants have anomalous distributions indicating that they are likely to be artefactual, `ALF` if relevant reads have lower median alignment score per base than a specified threshold, and `DVF` if variants appear to be the result of PCR error.

### FILTERS

- `ADF`; The `ADF` filter is an implementation of the artifact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads.

- `ALF`; The `ALF` filter indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` filter – artefacts with anomalous distributions often cause a marked decrease in alignment score.

- `DVF`; The `DVF` filter is a naive but effective algorithm for detecting variants which are the result of PCR error - in regions of low complexity, short repeats and homopolymer tracts can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling. The `DVF` filter attempts to catch these variants by examining the regularity of the start and end coordinates of collections of supporting reads and their mates.

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
The recommended usage is to provide a config of filter parameters along with the VCF in question and the relavant alignments (.sam/.bam/.cram), like so:
```
hairpin2 -c myconfig.json variants.vcf aln.cram > output.vcf
```
A config of default parameters is provided in `example-configs/`. All config parameters are equivalently named to their command line overrides (see helptext below), except `-` is replaced by `_`.

full helptext:
```
Usage: hairpin2 [-h, --help] [OPTIONS] VCF ALIGNMENTS...

  read-aware artefactual variant flagging algorithms. Flag variants in VCF
  using statistics calculated from supporting reads found in ALIGNMENTS, and
  emit the flagged VCF to stdout.

Options:
  -v, --version                   Show the version and exit.
  -h, --help                      show help (-h for basic, -hh for extended
                                  including config override options)
  -c, --config FILEPATH           path to config JSON from which filter
                                  paramters will be loaded - can be overridden
                                  by extended arguments provided at runtime
  -o, --output-config FILEPATH    log filter paramaters from run as a config
                                  JSON file
  -m, --name-mapping S:SM S:SM...
                                  If sample names in VCF differ from SM tags
                                  in alignment files, provide a key here to
                                  map them. When multiple alignments are
                                  provided, accepts a space separated list of
                                  sample:SM pairs. When only a single
                                  alignment is provided, also accepts a comma
                                  separated string of one or more possible
                                  sample-of-interest names like TUMOR,TUMOUR
  -r, --cram-reference FILEPATH   path to FASTA format CRAM reference,
                                  overrides $REF_PATH and UR tags for CRAM
                                  alignments
  -q, --quiet                     be quiet (-q to not log INFO level messages,
                                  -qq to additionally not log WARN)
  -p, --progress                  display progress bar on stderr during run

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
                                  inclusive maximum window size in number of
                                  bases within which read pairs supporting a
                                  variant may be considered possible
                                  duplicates of eachother. -1 will disable
                                  duplicate detection  [default: (6); x>=0]
    --loss-ratio FLOAT            ratio of the number of reads found to be
                                  duplicates against the total number of
                                  supporting reads, above which a variant is
                                  flagged DVF. In logical AND with a hardcoded
                                  test that at least 2 supporting reads are
                                  independent, i.e. not duplicates of each
                                  other, to ensure that regardless of the
                                  value of `loss_ratio` collapse of duplicates
                                  to only a single supporting read always
                                  results in a DVF flag. Smaller is more
                                  sensitive. Set to 0.99 to rely only on the
                                  hardcoded test (practically speaking).
                                  [default: (0.49); 0.0<=x<=0.99]

ALF config overrides:
    --al-filter-threshold FLOAT   threshold for median of read alignment score
                                  per base of all relevant reads, at and below
                                  which a variant is flagged ALF  [default:
                                  (0.93); 0<=x<=93]

ADF config overrides:
    --edge-definition FLOAT       percentage of a read that is considered to
                                  be "the edge" for the purposes of assessing
                                  variant position distribution  [default:
                                  (0.15); 0.0<=x<=0.99]
    --edge-fraction FLOAT         percentage of variants must occur within
                                  EDGE_FRACTION of read edges to mark ADF flag
                                  [default: (0.9); 0.0<=x<=0.99]
    --min-mad-one-strand INT      Mean Average Devaition of distances between
                                  variant position and read start above which
                                  a variant cannot be considered anomalous -
                                  used when only one strand has sufficient
                                  valid reads for testing  [default: (0);
                                  x>=0]
    --min-sd-one-strand FLOAT     stdev of distances between variant position
                                  and read start above which a variant cannot
                                  be considered anomalous - used when only one
                                  strand has sufficient valid reads for
                                  testing  [default: (4.0); x>=0.0]
    --min-mad-both-strand-weak INT
                                  Mean Average Devaition of distances between
                                  variant position and read start above which
                                  a variant cannot be considered anomalous -
                                  used when both strands have sufficient valid
                                  reads for testing, in logical AND with
                                  `min_sd_both_strand_weak`, and logical OR
                                  with corresponding strong condtion pair
                                  [default: (2); x>=0]
    --min-sd-both-strand-weak FLOAT
                                  stdev of distances between variant position
                                  and read start above which a variant cannot
                                  be considered anomalous - used when both
                                  strands have sufficient valid reads for
                                  testing, in logical AND with
                                  `min_mad_both_strand_weak`, and logical OR
                                  with corresponding strong condtion pair
                                  [default: (2.0); x>=0.0]
    --min-mad-both-strand-strong INT
                                  Mean Average Devaition of distances between
                                  variant position and read start above which
                                  a variant cannot be considered anomalous -
                                  used when both strands have sufficient valid
                                  reads for testing, in logical AND with
                                  `min_sd_both_strand_strong`, and logical OR
                                  with corresponding weak condtion pair
                                  [default: (1); x>=0]
    --min-sd-both-strand-strong FLOAT
                                  stdev of distances between variant position
                                  and read start above which a variant cannot
                                  be considered anomalous - used when both
                                  strands have sufficient valid reads for
                                  testing, in logical AND with
                                  `min_mad_both_strand_weak`, and logical OR
                                  with the corresponding weak condtion pair
                                  [default: (10.0); x>=0.0]
    --min-reads INT               number of reads at and below which the
                                  hairpin filtering logic considers a strand
                                  to have insufficient reads for testing for a
                                  given variant  [default: (1); x>=1]
```

Parameters are hopefully mostly clear from the helptext, but some warrant additional commentary:

- `--name-mapping` – When using multisample VCFS, hairpin2 compares VCF sample names found in the VCF header to SM tags in alignments to match samples of interest to the correct alignment. If these IDs are different between the VCF and alignments, you'll need to provide a key. If there are multiple samples of interest in the VCF, and therefore multiple alignments, you will need to provide a key for each pair - e.g. `-m sample1:SM1 sample2:SM2 ...`. If there is only one alignment, then you need only indicate which VCF sample is the sample of interest, e.g. `-m TUMOR`. As a convenience for high throughput workflows, when there is only one alignment you may also provide a comma separated string of possible names for the sample of interest, e.g. `-m TUMOR,TUMOUR`. Assuming there is one and only one match in the VCF, the tool will match the alignment to that sample.
- `--al-filter-threshold` – In the predecessor to `hairpin2`, `additionalBAMStatistics`, this value was known as `ASRD` and the default was set at 0.87.
- For all parameters, defaults were found by trial and error on LCM data and you may find it necessary to experiment with this parameter depending on data type.

Reading the `test()` method implementation of each filter may be informative. These can be found in `hairpin2/filters/`

The tool tests records in a VCF file and applies filter flags as appropriate. Reasoning for decisions is recorded in the INFO field of the VCF records, in the form `<FILTER>=<alt>|<True/False>|<code>|...`. Additional data (noted by the ellipsis) where present are described in the relevant header line of the VCF. The codes indicate the reason on which the decision was made, and are as follows:

DVF:
> **0** – insufficient supporting reads  
> **1** – variant is/not the result of PCR stuttering 

ALF:
> **0** – insufficient supporting reads  
> **1** – variant passed/failed on filter threshold 

ADF:
> **0** – insufficient supporting reads  
> **1** – passed/failed on condition 60A(i) of Ellis et al. (`ADF` only)  
> **2** – passed/failed on condition 60B(i) of Ellis et al. (`ADF` only)  

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
- stricter config and params validation, most likely with pydantic, to catch misformatted configs earlier
- further boundary condition testing
- improve documentation - describe filters in individual sections, beyond replicating helptext
- switch entirely to fstrings from .format()
- disscussions to be had on multisample VCF support


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
