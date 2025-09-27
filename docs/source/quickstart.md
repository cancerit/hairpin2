# Quickstart

### Installation
The easiest end-user approach is to install into a virtual environment:
```
python -m venv .env
source .env/bin/activate
pip install .
hairpin -h  # see helptext
```

for development, see README.md


### Usage

The recommended usage is to provide a config of flag parameters along with the VCF in question and the relavant alignments (.sam/.bam/.cram), like so:
```
hairpin2 -c myconfig.toml variants.vcf aln.cram > output.vcf
```
A config of default parameters is provided in `example-configs/default-params.toml`. Both TOML and JSON based configs are supported

