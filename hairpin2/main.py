# hairpin2
#
# Copyright (C) 2024, 2025 Genome Research Ltd.
#
# Author: Alex Byrne <ab63@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
from pathlib import Path
import pysam
from hairpin2 import  __version__
from hairpin2.abstractions.readawareproc import FlagResult
from hairpin2.const import TagEnum
from hairpin2.flaggers.LQF import FlaggerLQF, ResultLQF, TaggerLowQual
from hairpin2.flaggers.shared import RunParamsShared
from hairpin2.read_preprocessors.mark_overlap import TaggerOverlap
from hairpin2.read_preprocessors.mark_support import TaggerSupporting
from hairpin2.abstractions.structures import ReadView
from hairpin2.flaggers import ADF, ALF, DVF
import logging
import json
from itertools import tee
from typing import Any
import sys
import click
from dataclasses import dataclass, fields
# pyright: reportAny=false
# pyright: reportExplicitAny=false
# pyright: reportImplicitStringConcatenation=false

EXIT_SUCCESS = 0
EXIT_FAILURE = 1


@dataclass(frozen=True)
class DataclassTupleMixin:
    def __iter__(self):
        for f in fields(self):
            yield getattr(self, f.name)

    def __getitem__(self, index: int):
        return tuple(getattr(self, f.name) for f in fields(self))[index]


@dataclass(slots=True, frozen=True)
class MinMax(DataclassTupleMixin):
    min: int | float
    max: int | float | None = None


@dataclass(slots=True, frozen=True)
class ParamConstraint:
    default: int | float
    range: MinMax


_PARAMS = {
    'al_filter_threshold': ParamConstraint(0.93, MinMax(0,93)),
    'min_clip_quality': ParamConstraint(35, MinMax(0, 60)),
    'min_mapping_quality': ParamConstraint(11, MinMax(0, 93)),
    'min_base_quality': ParamConstraint(25, MinMax(-1)),
    'duplication_window_size': ParamConstraint(6, MinMax(0)),
    'loss_ratio': ParamConstraint(0.49, MinMax(0.0, 0.99)),
    'edge_definition': ParamConstraint(0.15, MinMax(0.0, 0.99)),
    'edge_fraction': ParamConstraint(0.9, MinMax(0.0, 0.99)),
    'min_mad_one_strand': ParamConstraint(0, MinMax(0)),
    'min_sd_one_strand': ParamConstraint(4.0, MinMax(0.0)),
    'min_mad_both_strand_weak': ParamConstraint(2, MinMax(0)),
    'min_sd_both_strand_weak': ParamConstraint(2.0, MinMax(0.0)),
    'min_mad_both_strand_strong': ParamConstraint(1, MinMax(0)),
    'min_sd_both_strand_strong': ParamConstraint(10.0, MinMax(0.0)),
    'min_reads': ParamConstraint(1, MinMax(1)),
}


existing_file_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
    resolve_path=True,
)

writeable_file_path = click.Path(
    dir_okay=False,
    writable=True,
    resolve_path=True,
)

class JSONOrFile(click.ParamType):
    def convert(self, value, param, ctx):
        # ty to interpret as path
        path_type = click.Path(exists=True, readable=True, dir_okay=False)
        try:
            file_path = path_type.convert(value, param, ctx)
            with open(file_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except click.BadParameter:
            pass  # Not a valid file path, treat as raw JSON
        except Exception as e:
            self.fail(f"Failed to parse JSON file '{value}': {e}", param, ctx)

        # Fallback: try parsing as raw JSON string
        try:
            return json.loads(value)
        except json.JSONDecodeError as e:
            self.fail(f"Invalid JSON input: {e}", param, ctx)



def show_help(ctx, value):
    if value > 1:
        # unhide hidden
        for param in ctx.command.params:
            if hasattr(param, "hidden"):
                param.hidden = False
            # for click-option-group: also unhide the parent group if needed
            if hasattr(param, "opt_group") and hasattr(param.opt_group, "hidden"):
                param.opt_group.hidden = False
    click.echo(ctx.get_help())
    ctx.exit()


@click.command(
    epilog='see documentation at https://github.com/cancerit/hairpin2 or at tool install location for further information',
    options_metavar='[-h, --help] [OPTIONS]',
    add_help_option=False,
)
@click.version_option(__version__, '-v', '--version', message='%(version)s')
@click.option(
    '-h',
    '--help',
    count=True,
    is_eager=True,
    expose_value=False,
    callback=lambda ctx, param, value: show_help(ctx, value) if value else None,
    help='show help (-h for basic, -hh for extended including config override options)'
)
@click.argument(
    'vcf',
    type=existing_file_path,
    required=True
)
@click.argument(
    'alignments',
    nargs=-1,
    type=existing_file_path,
    required=True
)

# Procedural options
@click.option(
    '-c',
    '--config',
    'config_path',
    metavar='FILEPATH',
    type=existing_file_path,
    help='path to config JSON from which filter paramters will be loaded - can be overridden by extended arguments provided at runtime'
)
@click.option(
    '-o',
    '--output-config',
    'output_config_path',
    metavar='FILEPATH',
    type=writeable_file_path,
    help='log filter paramaters from run as a config JSON file'
)
@click.option(
    '-m',
    '--name-mapping',
    metavar= "STR | FILEPATH",
    type=JSONOrFile(),
    help="If sample names in VCF differ from SM tags in alignment files, provide a key here to map them. "
    "Accepts a path to a JSON file, or JSON-formatted string of key-value pairs where keys are sample names in the VCF "
    "and all values are either the SM tag or the filepath of the relevant alignment "
    '- e.g. \'{"sample0": "PDxxA", "sample1": "PDxxB"}\' or \'{"sample0": "A.bam", ...}\'. '
    "When only a single alignment is provided, also accepts a JSON-spec top-level array of possible sample of interest names "
    '- e.g. \'["TUMOR","TUMOUR"]\'. '
    "Note that when providing a JSON-formatted string at the command line you must single quote the string, and use only double quotes internally"
)
@click.option(
    '-r',
    '--cram-reference',
    'cram_reference_path',
    metavar='FILEPATH',
    type=existing_file_path,
    help="path to FASTA format CRAM reference, overrides $REF_PATH and UR tags for CRAM alignments"
)
@click.option(
    '-q',
    '--quiet',
    count=True,
    help='be quiet (-q to not log INFO level messages, -qq to additionally not log WARN)',
    default=False
)
@click.option(
    '-p',
    '--progress',
    'progress_bar',
    is_flag=True,
    # metavar='',
    help='display progress bar on stderr during run',
    default=False
)

# # === Config override groups (hidden by default) ===
# @optgroup.group("\nread validation config overrides", hidden=True)
# @optgroup.option(
#     '--min-clip-quality',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_clip_quality'].range),
#     show_default=str(_PARAMS['min_clip_quality'].default),
#     help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases'
# )
# @optgroup.option(
#     '--min-mapping-quality',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_mapping_quality'].range),
#     show_default=str(_PARAMS['min_mapping_quality'].default),
#     help='discard reads with mapping quality below this value'
# )
# @optgroup.option(
#     '--min-base-quality',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_base_quality'].range),
#     show_default=str(_PARAMS['min_base_quality'].default),
#     help='discard reads with base quality below this value at variant position'
# )

# @optgroup.group('\nDVF config overrides', hidden=True)
# @optgroup.option(
#     '--duplication-window-size',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['duplication_window_size'].range),
#     show_default=str(_PARAMS['duplication_window_size'].default),
#     help='inclusive maximum window size in number of bases within which read pairs supporting a variant may be considered possible duplicates of eachother. -1 will disable duplicate detection'
# )
# @optgroup.option(
#     '--loss-ratio',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['loss_ratio'].range),
#     show_default=str(_PARAMS['loss_ratio'].default),
#     help='ratio of the number of reads found to be duplicates against the total number of supporting reads, above which a variant is flagged DVF. In logical AND with a hardcoded test that at least 2 supporting reads are independent, i.e. not duplicates of each other, to ensure that regardless of the value of `loss_ratio` collapse of duplicates to only a single supporting read always results in a DVF flag. Smaller is more sensitive. Set to 0.99 to rely only on the hardcoded test (practically speaking).'
# )

# @optgroup.group("\nALF config overrides", hidden=True)
# @optgroup.option(
#     '--al-filter-threshold',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['al_filter_threshold'].range),
#     show_default=str(_PARAMS['al_filter_threshold'].default),
#     help='threshold for median of read alignment score per base of all relevant reads, at and below which a variant is flagged ALF'
# )

# @optgroup.group("\nADF config overrides", hidden=True)
# @optgroup.option(
#     '--edge-definition',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['edge_definition'].range),
#     show_default=str(_PARAMS['edge_definition'].default),
#     help='percentage of a read that is considered to be "the edge" for the purposes of assessing variant position distribution'
# )
# @optgroup.option(
#     '--edge-fraction',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['edge_fraction'].range),
#     show_default=str(_PARAMS['edge_fraction'].default),
#     help='percentage of variants must occur within EDGE_FRACTION of read edges to mark ADF flag'
# )
# @optgroup.option(
#     '--min-mad-one-strand',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_mad_one_strand'].range),
#     show_default=str(_PARAMS['min_mad_one_strand'].default),
#     help='Mean Average Devaition of distances between variant position and read start above which a variant cannot be considered anomalous - used when only one strand has sufficient valid reads for testing'
# )  # 'above which' meaning exclusive
# @optgroup.option(
#     '--min-sd-one-strand',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['min_sd_one_strand'].range),
#     show_default=str(_PARAMS['min_sd_one_strand'].default),
#     help='stdev of distances between variant position and read start above which a variant cannot be considered anomalous - used when only one strand has sufficient valid reads for testing'
# )
# @optgroup.option(
#     '--min-mad-both-strand-weak',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_mad_both_strand_weak'].range),
#     show_default=str(_PARAMS['min_mad_both_strand_weak'].default),
#     help='Mean Average Devaition of distances between variant position and read start above which a variant cannot be considered anomalous - used when both strands have sufficient valid reads for testing, in logical AND with `min_sd_both_strand_weak`, and logical OR with corresponding strong condtion pair'
# )
# @optgroup.option(
#     '--min-sd-both-strand-weak',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['min_sd_both_strand_weak'].range),
#     show_default=str(_PARAMS['min_sd_both_strand_weak'].default),
#     help='stdev of distances between variant position and read start above which a variant cannot be considered anomalous - used when both strands have sufficient valid reads for testing, in logical AND with `min_mad_both_strand_weak`, and logical OR with corresponding strong condtion pair'
# )
# @optgroup.option(
#     '--min-mad-both-strand-strong',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_mad_both_strand_strong'].range),
#     show_default=str(_PARAMS['min_mad_both_strand_strong'].default),
#     help='Mean Average Devaition of distances between variant position and read start above which a variant cannot be considered anomalous - used when both strands have sufficient valid reads for testing, in logical AND with `min_sd_both_strand_strong`, and logical OR with corresponding weak condtion pair'
# )
# @optgroup.option(
#     '--min-sd-both-strand-strong',
#     metavar='FLOAT',
#     type=click.FloatRange(*_PARAMS['min_sd_both_strand_strong'].range),
#     show_default=str(_PARAMS['min_sd_both_strand_strong'].default),
#     help='stdev of distances between variant position and read start above which a variant cannot be considered anomalous - used when both strands have sufficient valid reads for testing, in logical AND with `min_mad_both_strand_weak`, and logical OR with the corresponding weak condtion pair'
# )
# @optgroup.option(
#     '--min-reads',
#     metavar='INT',
#     type=click.IntRange(*_PARAMS['min_reads'].range),
#     show_default=str(_PARAMS['min_reads'].default),
#     help='number of reads at and below which the hairpin filtering logic considers a strand to have insufficient reads for testing for a given variant'
# )
def hairpin2(
    vcf: str,
    alignments: list[str],
    config_path: str | None,
    output_config_path: str | None,
    name_mapping: dict[str, str] | list[str] | None,
    cram_reference_path: str | None,
    quiet: int,
    progress_bar: bool,
    **kwargs: Any,
) -> None:
    '''
    read-aware artefactual variant flagging algorithms. Flag variants in VCF using statistics calculated from supporting reads found in ALIGNMENTS, and emit the flagged VCF to stdout.
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s',
        datefmt='%I:%M:%S'
    )

    configd: dict[str, Any] = {}
    if config_path:
        try:
            with open(config_path, 'r') as f:
                configd= json.load(f)
        except Exception as er:
            logging.error(f'failed to open input JSON, reporting: {er}')
            sys.exit(EXIT_FAILURE)

    # verify config
    # for key in configd:
    #     if key not in _PARAMS:
    #         logging.error(f'unrecognised parameter {key!r} in config - check spelling and underscores?')
    #         sys.exit(EXIT_FAILURE)

    # override from config and ensure all args are set
    # test args are sensible, exit if not
    # unfortunately necessary duplication with options since the config remains unverified
    # TODO: move to validation of FilterParams with pydantic, where there should be full validation anyway
    # for key in kwargs:
    #     if kwargs[key] is None:
    #         if key in configd.keys():
    #             kwargs[key] = configd[key]
    #         elif key in _PARAMS:
    #             if not quiet: logging.info(f'parameter {key!r} not found in config or overridden on command line, falling back to standard default {_PARAMS[key].default}')
    #             kwargs[key] = _PARAMS[key].default
    #     if key not in _PARAMS:  # shouldn't happen
    #         logging.error(f'unrecognised parameter {key!r} in kwargs - this is probably a bug!')
    #         sys.exit(EXIT_FAILURE)
    #     else:
    #         assert kwargs[key] is not None
    #         if not isinstance(kwargs[key], type(_PARAMS[key].default)):
    #             logging.error(f'value {kwargs[key]!r} for parameter {key!r} is not of expected type {type(_PARAMS[key].default).__name__}')
    #             sys.exit(EXIT_FAILURE)
    #         pmin, pmax = _PARAMS[key].range
    #         if (pmin is not None and pmin > kwargs[key]) or (pmax is not None and kwargs[key] > pmax):
    #             logging.error(f'value {kwargs[key]!r} for parameter {key!r} out of range {pmin}<=x<={pmax}')
    #             sys.exit(EXIT_FAILURE)

    # # set up params based on inputs
    # shared_prefilter_d = {
    #     'min_mapq': kwargs['min_mapping_quality'],
    #     'min_avg_clipq': kwargs['min_clip_quality'],
    #     'min_baseq': kwargs['min_base_quality']
    # }

    # kwargs.update(shared_prefilter_d,)

    try:
        vcf_in_handle = pysam.VariantFile(vcf)
    except Exception as er:
        logging.error(f'failed to open input VCF, reporting {er}')
        sys.exit(EXIT_FAILURE)
    vcf_names: list[str] = list(vcf_in_handle.header.samples)
    if len(set(vcf_names)) != len(vcf_names):
        logging.error('duplicate sample names in input VCF')
        sys.exit(EXIT_FAILURE)

    sm_to_aln_map: dict[str, pysam.AlignmentFile] = {}
    filename_to_aln_map: dict[str, pysam.AlignmentFile] = {}
    for path_str in alignments:
        path = Path(path_str)
        try:
            match path.suffix[1]:
                case 's' | 'S':
                    mode = "r"
                case 'b' | 'B':
                    mode = "rb"
                case 'c' | 'C':
                    mode = "rc"
                case _:
                    raise ValueError
        except (IndexError, ValueError):
            logging.error(f'Could not infer alignment format from suffix {path.suffix!r} of file {path_str!r}')
            sys.exit(EXIT_FAILURE)
        if cram_reference_path and mode != 'rc':
            logging.error(f'CRAM reference provided, but alignment at {path_str!r} not inferred as cram from suffix {path.suffix!r}')
            sys.exit(EXIT_FAILURE)
        try:
            alignment = pysam.AlignmentFile(
                str(path),
                         mode,
                         reference_filename=(cram_reference_path if cram_reference_path else None)
            )
        except Exception as er:
            logging.error(f'failed to read alignment file {path!r}, reporting {er}')
            sys.exit(EXIT_FAILURE)
        # grab the sample name from first SM field
        # in header field RG
        aln_sm = alignment.header.to_dict()['RG'][0]['SM']
        sm_to_aln_map[aln_sm] = alignment
        filename_to_aln_map[path.name] = alignment

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}

    match name_mapping:
        case list():
            matches = [name for name in name_mapping if name in set(vcf_names)]
            if len(matches) > 1:
                logging.error(msg='More than one of the VCF sample names provided to name mapping {name_mapping} match any sample names in input VCF {vcf_names}!')
                sys.exit(EXIT_FAILURE)
            elif not matches:
                logging.error(msg=f'None of VCF sample names provided to name mapping {name_mapping} match any sample name in input VCF {vcf_names}!')
                sys.exit(EXIT_FAILURE)
            else:
                if not quiet: logging.info('matched alignment to sample {} in VCF'.format(matches[0]))
                vcf_sample_to_alignment_map[matches[0]] = alignment  # pyright: ignore[reportPossiblyUnboundVariable] | since length of alignments == 1, can just reuse this variable
        case dict():
            vcf_map_keys: list[str] = []
            alignment_map_values: list[str] = []
            for key, val in name_mapping.items():
                vcf_map_keys.append(key)
                alignment_map_values.append(val)

            if len(vcf_map_keys) != len(set(vcf_map_keys)):
                logging.error('duplicate keys (VCF sample names) provided to name mapping!')
                sys.exit(EXIT_FAILURE)

            if not set(vcf_map_keys) <= set(vcf_names):
                logging.error(f'keys (VCF sample names) provided to name mapping {vcf_map_keys!r} are not equal to or a subset of VCF samples from input VCF {vcf_names!r}!')
                sys.exit(EXIT_FAILURE)

            if len(alignment_map_values) != len(set(alignment_map_values)):
                logging.error(msg='duplicate values (alignment SM tags/filenames) provided to name mapping!')
                sys.exit(EXIT_FAILURE)

            match_sm = sorted(alignment_map_values) == sorted(sm_to_aln_map.keys())
            match_fn = sorted(alignment_map_values) == sorted(filename_to_aln_map.keys())
            if match_sm:
                vcf_sample_to_alignment_map = {
                    vcf_map_keys[alignment_map_values.index(k)]: v
                    for k, v
                    in sm_to_aln_map.items()
                }
            elif match_fn:
                vcf_sample_to_alignment_map = {
                    vcf_map_keys[alignment_map_values.index(k)]: v
                    for k, v
                    in filename_to_aln_map.items()
                }
            else:
                logging.error(msg=f'values provided to name mapping {alignment_map_values!r} are not equal to the either the SM tags or the filenames of the alignments provided')
                sys.exit(EXIT_FAILURE)

            if match_sm and match_fn:
                logging.warning(f'intention ambigous - values provided to name mapping {alignment_map_values!r} are equal to both the SM tags and the filenames of alignments. Mapping against SM tags, not filenames')
        case None:
            if not sm_to_aln_map.keys() <= set(vcf_names):
                logging.error(
                    msg='alignment SM tags {} are not equal to or a subset of VCF sample names {} - use a name mapping'.format(
                        set(sm_to_aln_map.keys()),
                        set(vcf_names)
                    )
                )
                sys.exit(EXIT_FAILURE)
            vcf_sample_to_alignment_map = sm_to_aln_map
        case _:
            logging.error(f'name mapping misformatted. Expected deserialised JSON, recieved {type(name_mapping)}, containing {name_mapping}')
            sys.exit(EXIT_FAILURE)

    if set(vcf_names) != vcf_sample_to_alignment_map.keys():
        if not quiet: logging.info(
            "alignments not provided for all VCF samples; {} will be ignored".format(
                set(vcf_names) - vcf_sample_to_alignment_map.keys()
            )
        )

    # init output  TODO: put these in filter modules as funcs
    out_head = vcf_in_handle.header
    out_head.add_line(
        f'##FILTER=<ID=ALF,Description="Median alignment score of reads reporting variant less than {configd['ALF']['avg_AS_threshold']}">'  # TODO: move to result or flagger classes
    )
    out_head.add_line(
        f'##FILTER=<ID=ADF,Description="Variant shows anomalous distribution in supporting reads">'
    )
    out_head.add_line(
        f'##FILTER=<ID=DVF,Description="More than {configd['DVF']['read_loss_threshold']} of reads supporting this variant are considered PCR stutter duplicates">'
    )
    out_head.add_line(
        f'##FILTER=<ID=LQF,Description="More than {configd['LQF']['read_loss_threshold']} of reads supporting this variant are considered low quality">'
    )
    out_head.add_line(
        '##INFO=<ID=ADF,Number=1,Type=String,Description="alt|[True,False]|code indicating decision reason for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=ALF,Number=1,Type=String,Description="alt|[True,False]|code|score indicating decision reason and average AS of supporting reads for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=DVF,Number=1,Type=String,Description="alt|[True,False]|code|loss indicating decision reason for each alt and ratio of supporting reads suspected to be duplicates">'
    )
    out_head.add_line(
        '##INFO=<ID=LQF,Number=1,Type=String,Description="alt|[True,False]|code|loss indicating decision reason for each alt and ratio of supporting reads suspected to be low quality">'
    )
    out_head.add_line(
        f'##hairpin2_version={__version__}'
    )
    out_head.add_line(
        f'##hairpin2_params=[{json.dumps(configd)}]'
    )
    out_head.add_line(
        f'##hairpin2_samples={vcf_sample_to_alignment_map.keys()}'
    )

    try:
        vcf_out_handle = pysam.VariantFile(sys.stdout, 'w', header=out_head)
    except Exception as e:
        logging.error(msg='failed to open VCF output, reporting: {}'.format(e))
        sys.exit(EXIT_FAILURE)


    # NOTE: all this below should be excised into raf eventually
    sp = TaggerSupporting(configd, [], [])
    ov = TaggerOverlap(configd, [TagEnum.SUPPORT], [])
    lqt = TaggerLowQual(configd, [TagEnum.SUPPORT], [])
    dp = DVF.TaggerDupmark(configd, [TagEnum.SUPPORT], [TagEnum.LOW_QUAL])

    lqf = FlaggerLQF(configd, [TagEnum.SUPPORT], [TagEnum.OVERLAP])
    ad = ADF.FlaggerADF(configd, [TagEnum.SUPPORT], [TagEnum.STUTTER_DUP, TagEnum.OVERLAP, TagEnum.LOW_QUAL])
    al = ALF.FlaggerALF(configd, [TagEnum.SUPPORT], [TagEnum.OVERLAP, TagEnum.LOW_QUAL, TagEnum.STUTTER_DUP])
    dv = DVF.FlaggerDVF(configd, [TagEnum.SUPPORT], [TagEnum.OVERLAP, TagEnum.LOW_QUAL])
    # test records
    prog_bar_counter = 0
    for record in vcf_in_handle.fetch():
        record_filters: dict[type[FlagResult[Any]], list[FlagResult[Any]]] = {ALF.ResultALF: [], ADF.ResultADF: [], DVF.ResultDVF: [], ResultLQF: []}
        if record.alts is None:
            if quiet < 1: logging.warning('{0: <7}:{1: >12} ¦ no alts for this record, skipping'.format(
                record.chrom, record.pos))
        else:
            # TODO: also need to mandate/require field presence on variant records I suppose
            samples_w_mutants = [name
                                 for name
                                 in record.samples
                                 if record.samples[name]["GT"] != (0, 0)]
            if len(samples_w_mutants) == 0:
                if quiet < 1: logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit alts associated with this record, skipping'.format(
                    record.chrom, record.pos))
            else:
                reads_by_sample: dict[str, list[pysam.AlignedSegment]] = {}

                # get pileup
                for k, v in vcf_sample_to_alignment_map.items():
                    if k in samples_w_mutants:
                        read_iter, test_iter = tee(v.fetch(record.chrom,
                                                           record.start,
                                                           (record.start + 1)))
                        try:
                            _ = next(test_iter)
                        except StopIteration:
                            continue  # no reads for that sample cover this region
                        else:
                            reads_by_sample[k] = list(read_iter)

                # --- test by alt ---
                # TODO: put mutation type detection under testing
                # the ability to handle complex mutations would be a potentially interesting future feature
                # for extending to more varied artifacts
                # TODO/NOTE: this should itself be an additive process! on variant rather than reads
                for alt in record.alts:
                    if (record.rlen == len(alt)
                            and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):
                        mut_type = 'S'
                    elif (len(alt) < record.rlen
                            and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):  # DEL - DOES NOT SUPPORT <DEL> TYPE IDS OR .
                        mut_type = 'D'
                    elif (record.rlen == 1
                            and set(alt).issubset(set(['A', 'C', 'T', 'G', 'N', '*']))):  # INS - DOES NOT SUPPORT <INS> TYPE IDS
                        mut_type = 'I'
                    else:
                        if quiet < 1: logging.warning('could not infer mutation type, POS={} REF={} ALT={}, skipping variant'.format(
                            record.pos, record.ref, alt))
                        continue

                    # FUTURE: shared qc filtering based on parsing of flagger prefilter configs, post additive preprocessing
                    # NOTE: global check for tag availability on reads, and filtering of reads without those properties
                    # would save a lot of performance as fields checked globally could be subtracted from fields checked
                    # on specific flaggers and preprocessors, and less reads to hold.


                    # instantiate test data obj/s
                    # NOTE: shared var params KOs flagger independence
                    # NOTE: shared access to the same record in memory is dangerous, trusts developers not to modify for all
                    # in flagger methods
                    # NOTE: doubly so if multithreading in future
                    # reasonably immutable view of test data
                    # if further immutablility is needed will have to wrap alignedsegment internally when priming
                    # TODO/NOTE: wrapping pysam types with wrapt or similar for more defined behaviour with method forwarding
                    # and tag registries etc
                    # TODO/NOTE: dict-by-tag bevhaviour in ReadView will allow for really neat subselection of reads based on
                    # intersection of include/exclude tags
                    test_reads = ReadView(ReadView.convert_pysam_to_extread(reads_by_sample, validate=True))
                    run_data = RunParamsShared(record=record, reads=test_reads, alt=alt, mut_type=mut_type)  # TODO: allow positional args
                    # dpp(dupmark.RunParamsDupmark(record=record, reads=test_reads, alt=alt, mut_type=mut_type))

                    # run the flaggers
                    # TODO: run all via centralised runner


                    # SUCCESS! NO REGRESSIONS AGAINST 2.0.0 ON MY DATA. TODO: separate args, wire up to CLI interface properly
                    sp(run_data)
                    ov(run_data)
                    lqt(run_data)
                    dp(run_data)  # needs lqt run first - TODO: this fact isn't surfaced or guarded in any way
                    # TODO: to test this, need to assess passing reads using these tag based methods vs passing reads in old qc funcs implementation. If they're the same
                    # then the only remaining problem space is the LQF implementation
                    lq_result = lqf(run_data)
                    dv_result = dv(run_data)
                    al_result = al(run_data)
                    ad_result = ad(run_data)
                    for res in (al_result, ad_result, dv_result, lq_result):
                        assert res is not None
                        record_filters[type(res)].append(res)

                    sp.reset()
                    ov.reset()
                    lqt.reset()
                    dp.reset()
                    lqf.reset()
                    dv.reset(); al.reset(); ad.reset()  # clear for next variant

        if any(lst for lst in record_filters.values()):
            for ftype in record_filters:
                if any(fres.flag == True for fres in record_filters[ftype]):
                    record.filter.add(ftype.FlagName)
                # TODO: LQF prints no info!
                record.info.update({ftype.FlagName: ','.join([fl.getinfo() for fl in record_filters[ftype]])})  # pyright: ignore[reportArgumentType]

        try:
            vcf_out_handle.write(record)
        except Exception as e:
            logging.error(msg='failed to write to vcf, reporting: {}'.format(e))
            sys.exit(EXIT_FAILURE)

        if progress_bar:
            if not prog_bar_counter % 100:
                print(f'\rchr: {record.chrom}, pos: {record.pos}', end='', flush=True, file=sys.stderr)
            prog_bar_counter += 1
    else:
        if progress_bar:
            print(file=sys.stderr)

    if output_config_path:
        try:
            with open(output_config_path, "w") as output_json:
                json.dump(
                    {
                        k: configd[k]
                        for k
                        in _PARAMS
                    },
                    output_json, indent="")
        except Exception as e:
            logging.error(msg='failed to write output JSON, reporting: {}'.format(e))
            sys.exit(EXIT_FAILURE)

    if not quiet: logging.info('hairpin complete')
    sys.exit(EXIT_SUCCESS)

