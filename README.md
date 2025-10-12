[![python](https://img.shields.io/badge/Python-3.12-blue?style=for-the-badge&logo=python&logoColor=FFD43B)](https://docs.python.org/3.12/)
[![license](https://img.shields.io/badge/License-MIT-a51931?style=for-the-badge)](LICENSE.txt)
[![Documentation](https://img.shields.io/badge/Documentation-Online-blue?style=for-the-badge&logo=readthedocs)](https://cancerit.github.io/hairpin2/)

# hairpin2

[DOCUMENTATION HERE](https://cancerit.github.io/hairpin2/)

`hairpin2` – read-aware artefactual variant flagging

`hairpin2` is designed to flag variants that are likely artefactual via a series of tests performed upon the read data associated with each variant. Initially, it was concieved to flag possible cruciform artefacts for LCM sequence data, but the concept has been extended and can detect a variety of potentially spurious variants (including indels). The tool operates on a VCF file containing one or more samples, and alignment files for all samples to be tested.

Given a VCF, and BAM files for the samples of that VCF, return a VCF with variants flagged with `ADF` if variants have anomalous distributions indicating that they are likely to be artefactual, `ALF` if relevant reads have lower median alignment score per base than a specified threshold, `DVF` if variants appear to be the result of PCR error, and `LQF` if the variant is largely supported by low quality reads.

### FLAGS

- `ADF`; The `ADF` flag is an implementation of the artifact detection algorithm described in [Ellis et al, 2020](https://www.nature.com/articles/s41596-020-00437-6). It detects variants which appear with anomalously regular positional distribution in supporting reads.

- `ALF`; The `ALF` flag indicates variants which are supported by reads with poor signal-to-noise, per the alignment score. It is complementary to the `ADF` flag – artefacts with anomalous distributions often cause a marked decrease in alignment score.

- `DVF`; The `DVF` flag is a naive but effective algorithm for detecting variants which are the result of PCR error - in regions of low complexity, short repeats and homopolymer tracts can cause PCR stuttering. PCR stuttering can lead to, for example, an erroneous additional A on the read when amplifying a tract of As. If duplicated reads contain stutter, this can lead to variation of read length and alignment to reference between reads that are in fact duplicates. Because of this, these duplicates both evade dupmarking and give rise to spurious variants when calling. The `DVF` flag attempts to catch these variants by examining the regularity of the start and end coordinates of collections of supporting reads and their mates.

 - `LQF`; The `LQF` flag is a superset of the `DVF` flag - it tests whether a read is largely supported by both low quality reads and stutter duplicate reads (which are also considered low quality). Note that because the parameters for each `LQF` and `DVF` are independent, you can indepedently set the sensitivity of each - so the result of LQF is not necessarily a complete overlap with DVF (and usually is not).  

All flags are tunable such that their parameters can be configured to a variety of use cases and sequencing methods.

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
pip install -e ".[dev,doc]"
```

### ASSUMPTIONS & LIMITATIONS
`hairpin2` is designed for paired data where alignment records have the `MC` tag and the complete CIGAR string is present in the `CIGAR` field (rather than the `CG:B,I` tag). If the `MC` tag is not present in your data, it can be added using `samtools fixmate` or `biobambam2 bamsormadup`. The tool can handle substitions, insertions, and deletions formatted per the VCF specification. At this time, the tool will not investigate mutations notated with angle brackets, e.g. `<DEL>`, complex mutations, or monomorphic reference. No further assumptions are made – other alignment tags and VCF fields are used, however they are mandatory per the relevant format specifications. If these requirements are limiting and you need the tool to be extended in some way, please request it.

<!-- drop until testing reinstated -->
<!-- ### TESTING -->
<!-- A test suite has been provided with the algorithm implementation. To run these tests run `pytest .` from within the install directory. `hairpin2` must have been installed from that same directory, and be available on path (for example in a virtual environment). The tests can be found in the `test` directory. The focus is on testing all nodes (statements) and edges (paths between statements) for key scientific logic, rather than trying every possible input combination, in addition to testing critical boundary conditions. Since all input possibilities will pass via the same network/graph, if we prove that each part of that network functions correctly, we can be more confident that the program functions as it claims to. The tests are simple, and, once you have read them, it should be very easy to add your own further tests should you feel the need to confirm any behaviour. -->


### TODO
- automated regression testing
- disscussions to be had on multisample VCF support, and multiallelic variant support

### AUTHORS

`hairpin2` was developed at the Wellcome Sanger Institue via collaboration between CASM Informatics and CASM faculty. We thank all contributors and scientific collaborators for their input and expertise.

`hairpin2` is chiefly the work of:

Alex Byrne - [blex.bio](blex.bio) - Lead Developer/Contact  
Anh Phuong Le - [GitHub](https://github.com/Phuong-Le)  
Peter Campbell - [Quotient Therapeutics](https://quotient-tx.com/team/peter-campbell)

