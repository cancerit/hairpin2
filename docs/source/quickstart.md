(quickstart-target)=
# Quickstart

### Installation
The easiest end-user approach is to install into a virtual environment:
```bash
python -m venv .env
source .env/bin/activate
pip install .
hairpin --help
```

### Usage

The recommended usage is to provide a config of flag parameters along with the VCF in question and the relavant alignment/s (.sam/.bam/.cram), like so:
```bash
vcf_sample_name="TUMOUR"
aln="aln.cram"
hairpin2 \
  -c myconfig.toml \
  -m '{"$vcf_sample_name":"$aln"}' \
  variants.vcf \
  $aln > output.vcf
```
A config of default parameters is provided in `example-configs/default-params.toml`. Both TOML and JSON based configs are supported.

The default parameters provided are those found to be appropriate on series of LCMB data. Your use case, data, and opinions may differ - hairpin2 is extensively customisable via the config, and descriptions of the parameters can be found in the [Processes and Parameters](#process-target) section of the guide.  

See the [guide](#guide-target) for a complete walkthrough, and the [interface](#Command_Line_Interface) section below for a description of the command line options.

### Command Line Interface

```
Usage: hairpin2 [-h, --help] [OPTIONS] VCF ALIGNMENTS...

  read-aware artefactual variant flagging algorithms. Flag variants in VCF
  using statistics calculated from supporting reads found in ALIGNMENTS and
  emit the flagged VCF to stdout.

Usage: hairpin2 [-h, --help] [OPTIONS] VCF ALIGNMENTS...

  read-aware artefactual variant flagging algorithms. Flag variants in VCF
  using statistics calculated from supporting reads found in ALIGNMENTS and
  emit the flagged VCF to stdout.

Options:
  -v, --version                   Show the version and exit.
  -c, --config FILEPATH           path to config TOML/s or JSON/s from which
                                  processes and execution will be configured.
                                  May be provided multiple times; individual
                                  configs can be provided for each top level
                                  key (params, exec).
  -o, --output-config FILEPATH    log run configuration back to a new JSON
                                  file.
  --set <TEXT TEXT>...            Override values in the supplied config. Uses
                                  dot paths for the key, and JSON format for
                                  the value, e.g., --set
                                  params.DVF.read_loss_threshold 0.6. Must be
                                  provided after --config. May be provided
                                  multiple times.
  -m, --name-mapping JSON_STRING | FILEPATH
                                  If sample names in VCF differ from SM tags
                                  in alignment files, provide a key here to
                                  map them. Accepts a path to a JSON file, or
                                  JSON-formatted string of key-value pairs
                                  where keys are sample names in the VCF and
                                  all values are either the SM tag or the
                                  filepath of the relevant alignment - e.g.
                                  '{"sample0": "PDxxA", "sample1": "PDxxB"}'
                                  or '{"sample0": "A.bam", ...}'. When only a
                                  single alignment is provided, also accepts a
                                  JSON-spec top-level array of possible sample
                                  of interest names - e.g.
                                  '["TUMOR","TUMOUR"]'. Note that when
                                  providing a JSON-formatted string at the
                                  command line you must single quote the
                                  string, and use only double quotes
                                  internally.
  -r, --cram-reference FILEPATH   path to FASTA format CRAM reference,
                                  overrides $REF_PATH and UR tags for CRAM
                                  alignments.
  -q, --quiet                     be quiet (-q to not log INFO level messages,
                                  -qq to additionally not log WARN).
  -p, --progress                  display progress bar on stderr during run.
  -h, --help                      Show this message and exit.
```

To expand on `--name-mapping` – when using multisample VCFs, hairpin2 compares VCF sample names found in the VCF header to SM tags in alignments to match samples of interest to the correct alignment. If these IDs are different between the VCF and alignments, you'll need to provide a JSON key. If there are multiple samples of interest in a multisample VCF, and therefore it is necessary to provide multiple alignments, you will need to provide a mapping for each pair - e.g. `-m '{"sample1":"SM1", "sample2":"SM2", ...}'` or `-m '{"sample1:"1.bam", ...}'`. If there is only one sample of interest, and therefore only one alignment is provided to the tool, then you also have an optional shorthand - you need only indicate which VCF sample is the sample of interest, e.g. `-m '["TUMOR"]'`. When there is only one sample of interest, and therefore one alignment, but the sample of interest may have one of several possible names, you may also provide a comma separated string of possible names for the sample of interest, e.g. `-m '["TUMOR", "TUMOUR"]'` - users have found this valuable for high throughput workflows where the VCF input may be coming from one of several prior tools or callers (which may name samples differently). In all cases, there must be one and only one match between each alignment and VCF sample. In all cases, a path to a JSON file may be provided instead of the JSON string. Note that a VCF containing both a TUMOUR and a NORMAL sample contains 2 samples, and therefore is a multisample VCF.


### hp2-utils

`hairpin2` comes with a daughter tool, `hp2-utils`. `hp2-utils` provides some additional functionality to assist with the usage of hairpin2.  

```
Usage: hp2-utils [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  explain-var  explain a hairpin2 flagging decision
  get-params   get run parameters from a VCF
```

(explain-qs-target)=
#### explain-var
`explain-var` operates on any of the INFO fields set by hairpin2 in the output VCF, revealing the specific reasoning as to the decision the tool made for a given flag. This is very useful when tuning hairpin2 to new data, or querying results about which you are curious. Usage is as follows:
```bash
hp2-utils explain-var "ADF=G|PASS|0x77|19|BOTH"  # an example INFO string
```
returning:
```
FLAG: ADF
ALT: G
variant outcome was PASS via conditions ['NO_TESTABLE_READS', 'INSUFFICIENT_READS', 'EDGE_CLUSTERING', 'BOTH_STRAND_DISTRIB_BOTH', 'BOTH_STRAND_DISTRIB_ONE', 'MIN_NON_EDGE'] on strand BOTH
reads examined: 19
```

More information can be found at [Understanding Decisions](#understanding-decisions-target).

(get-params-qs-target)=
#### get-params
hairpin2 stores all run parameters in a reduced representation in the header of the output VCF. `get-params` allows translation of the stored parameters back into a JSON, which can then be used for other hairpin2 runs. Usage of `get-params` is as follows:
```bash
hp2-utils get-params my.hairpin2.vcf > parameters_used.json
```

More information can be found at [Reproducibility & Parameter Distribution](#repro-target).
