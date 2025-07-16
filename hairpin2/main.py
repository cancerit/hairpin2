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
from hairpin2.abstractfilters import FilterResult
from hairpin2.filters import ADF, ALF, DVF
from hairpin2.readqc import qc_read_broad, qc_read_alt_specific
import logging
import json
from itertools import tee, chain
from typing import Any
from collections.abc import Iterable
import sys
import click
from click_option_group import optgroup
from enum import Enum
from dataclasses import dataclass, fields
# pyright: reportUnusedCallResult=false
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
    min: int
    max: int


class _DEFAULTS(Enum):
    al_filter_threshold = 0.93
    min_clip_quality = 35
    min_mapping_quality = 11
    min_base_quality = 25
    duplication_window_size = 6
    edge_definition = 0.15
    edge_fraction = 0.9
    min_MAD_one_strand = 0
    min_sd_one_strand = 4.0
    min_MAD_both_strand_weak = 2
    min_sd_both_strand_weak = 2.0
    min_MAD_both_strand_strong = 1
    min_sd_both_strand_strong = 10.0
    min_reads = 1


# TODO: parameter _RANGES for read qc params
# TODO: validation of appropriate params in FilterParams (use pydantic)
class _RANGES(Enum):
    ...


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


def show_full_help(ctx):
    # Unhide all hidden option groups and options
    for param in ctx.command.params:
        if hasattr(param, "hidden"):
            param.hidden = False
        # for click-option-group: also unhide the parent group if needed
        if hasattr(param, "opt_group") and hasattr(param.opt_group, "hidden"):
            param.opt_group.hidden = False

    click.echo(ctx.get_help())
    ctx.exit()


@click.command(
    no_args_is_help=True,
    epilog='see documentation at https://github.com/cancerit/hairpin2 or at tool install location for further information',
    options_metavar='[-h, --help] [OPTIONS]')
@click.version_option(__version__, '-v', '--version', message='%(version)s')
@click.help_option('-h', '--help', help='show helptext')
@click.option('--help-all', '-hh', is_flag=True, is_eager=False, expose_value=False,
              callback=lambda ctx, param, value: show_full_help(ctx),
              help='Show full help including config override options.')
@click.argument('vcf_in', type=existing_file_path, required=True)
@click.argument('alignments', nargs=-1, type=existing_file_path, required=True)
@click.option('-f', '--format', type=click.Choice(['s', 'b', 'c']), required=True,
              help="format of alignment files; s indicates SAM, b indicates BAM, and c indicates CRAM")  # TODO: surely we can infer this

# Procedural options
@click.option('-r', '--cram-reference', metavar='FILEPATH', type=existing_file_path,
              help="path to FASTA format CRAM reference, overrides $REF_PATH and UR tags - ignored if --format is not CRAM")
@click.option('-m', '--name-mapping', nargs=0, multiple=True,
              help="key to map samples in a multisample VCF to alignment/s provided to -a. Uses VCF sample names per VCF header "
                   "and alignment SM tags. With multiple alignments to -a, accepts a space separated list of sample:SM pairs. "
                   "With a single alignment, also accepts a comma separated string of one or more possible sample-of-interest "
                   "names like TUMOR,TUMOUR")  # TODO: metavar
@click.option('-c', '--config', 'config_path', metavar='FILEPATH', type=existing_file_path,  # Validate config and populate subsequent args with click (I think this is a thing)
              help='path to config JSON from which filter paramters will be loaded - can be overridden by extended arguments provided at runtime')
@click.option('-o' , '--output_config', 'output_config_path', metavar='FILEPATH', type=writeable_file_path,
              help='log filter paramaters from run as a config JSON file')
@click.option('--progess/--no-progess', 'progress_bar', help='display progress bar during run', default=False)

# === Config override groups (hidden by default) ===
@optgroup.group("\nread validation config overrides", hidden=True)
@optgroup.option('--min-clip-quality', metavar='INT', type=click.IntRange(0, 93), default=None, show_default=str(_DEFAULTS.min_clip_quality),
                 help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases')
@optgroup.option('--min-mapping-quality', metavar='INT', type=click.IntRange(0, 60), default=11, show_default=str(_DEFAULTS.min_mapping_quality),
                 help='discard reads with mapping quality below this value')
@optgroup.option('--min-base-quality', metavar='INT', type=click.IntRange(0, 93), default=25, show_default=str(_DEFAULTS.min_base_quality),
                 help='discard reads with base quality below this value at variant position')

@optgroup.group('DVF config overrides', hidden=True)
@optgroup.option('--duplication-window-size', metavar='INT', type=click.IntRange(-1), default=6, show_default=str(_DEFAULTS.duplication_window_size),
                 help='inclusive maximum window size, in number of bases to use when detecting PCR duplicates. -1 will disable duplicate detection')  # TODO/BUG: is -1 still true?

@optgroup.group("\nALF config overrides", hidden=True)
@optgroup.option('--al-filter-threshold', metavar='FLOAT', type=click.FloatRange(0.0), default=0.93, show_default=str(_DEFAULTS.al_filter_threshold),
                 help='ALF; threshold for median of read alignment score per base of all relevant reads, at and below which a variant is flagged as ALF')

@optgroup.group("\nADF config overrides", hidden=True)
@optgroup.option('--edge-definition', metavar='FLOAT', type=click.FloatRange(0.0, 1.0), default=0.15, show_default=str(_DEFAULTS.edge_definition),
                 help='ADF; percentage of a read that is considered to be "the edge" for the purposes of assessing variant location distribution')
@optgroup.option('--edge-fraction', metavar='FLOAT', type=click.FloatRange(0.0, 1.0), default=0.9, show_default=str(_DEFAULTS.edge_fraction),
                 help='ADF; percentage of variants must occur within EDGE_FRACTION of read edges to mark ADF flag')
@optgroup.option('--min-mad-one-strand', metavar='INT', type=click.IntRange(0), default=0, show_default=str(_DEFAULTS.min_MAD_one_strand),
                 help='ADF; min range of distances between variant position and read start for valid reads when only one strand has sufficient valid reads for testing')  # BUG: what?
@optgroup.option('--min-sd-one-strand', metavar='FLOAT', type=click.FloatRange(0.0), default=4.0, show_default=str(_DEFAULTS.min_sd_one_strand),
                 help='ADF; min stdev of variant position and read start for valid reads when only one strand has sufficient valid reads for testing') # BUG: exclusive?!
@optgroup.option('--min-mad-both-strand-weak', metavar='INT', type=click.IntRange(0), default=2, show_default=str(_DEFAULTS.min_MAD_both_strand_weak),
                 help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbsw is true')
@optgroup.option('--min-sd-both-strand-weak', metavar='FLOAT', type=click.FloatRange(0.0), default=2.0, show_default=str(_DEFAULTS.min_sd_both_strand_weak),
                 help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbsw is true- default: 2, range: 0-, exclusive')
@optgroup.option('--min-mad-both-strand-strong', metavar='INT', type=click.IntRange(0), default=1, show_default=str(_DEFAULTS.min_MAD_both_strand_strong),
                 help='ADF; min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbss is true')
@optgroup.option('--min-sd-both-strand-strong', metavar='FLOAT', type=click.FloatRange(0.0), default=10, show_default=str(_DEFAULTS.min_sd_both_strand_strong),  # TODO: need a defaults dataclass or constants
                 help='ADF; min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbss is true')
@optgroup.option('--min-reads', metavar='INT', type=click.IntRange(0), default=min, show_default=_DEFAULTS.min_reads,
                 help='ADF; number of reads at and below which the hairpin filtering logic considers a strand to have insufficient reads for testing')
# TODO: DVF option decorators!
def hairpin2(
    input_vcf_path: str,
    alignment_paths: list[str],
    format: str,
    config_path: str | None = None,  # TODO: take json, load via click decorator?
    progress_bar: bool = False,
    cram_reference_path: str | None = None,
    config_output_path: Path | None = None,
    **kwargs: Any,
    # name_mapping: list[str] | None = None,
    # min_clip_quality: int | None = None,
    # min_mapping_quality: int | None = None,
    # min_base_quality: int | None = None,
    # duplication_window_size: int | None = None,
    # al_filter_threshold:  float | None = None,
    # edge_definition: float | None = None,
    # edge_fraction: float | None = None,
    # min_mad_one_strand: int | None = None,
    # min_sd_one_strand: float | None = None,
    # min_mad_both_strand_weak: int | None = None,
    # min_sd_both_strand_weak: float | None = None,
    # min_mad_both_strand_strong: int | None = None,
    # min_sd_both_strand_strong: float | None = None,
    # min_reads: int | None = None,
) -> None:
    '''
    cruciform artefact flagging algorithm based on Ellis et al. 2020 (DOI: 10.1038/s41596-020-00437-6)
    '''
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s ¦ %(levelname)-8s ¦ %(message)s',
        datefmt='%I:%M:%S'
    ) # TODO: optionally log to file since this may contain warnings about variants? Or to VCF...?

    configd: dict[str, Any] = {}
    if config_path:
        try:
            with open(config_path, 'r') as f:
                configd= json.load(f)
        except Exception as er:
            logging.error(f'failed to open input JSON, reporting: {er}')
            sys.exit(EXIT_FAILURE)

    # override from config and ensure all args are set
    for key in kwargs:
        if not kwargs[key]:
            if key in configd.keys():
                setattr(kwargs, key, configd[key])
            elif str(key) in _DEFAULTS.__members__:
                logging.info(f'arg {key!r} not found in config or present on command line, falling back to standard default {_DEFAULTS(key)}')
                setattr(kwargs, key, _DEFAULTS(key))
            elif str(key) in _DEFAULTS.__members__:
                setattr(kwargs, key, _DEFAULTS(key))

    # test args are sensible, exit if not
    # unfortunately necessary duplication with options since the config remains unverified
    # but can move most of these to validation of FilterParams, where there should be full validation anyway
    if not any([(0 <= min_clip_quality <= 93),
                (0 <= min_mapping_quality <= 60),
                (0 <= min_base_quality <= 93),
                (duplication_window_size >= -1),
                (al_filter_threshold >= 0),
                (0 <= edge_definition <= 0.99),
                (0 <= edge_fraction <= 0.99),
                (min_mad_one_strand >= 0),
                (min_sd_one_strand >= 0),
                (min_mad_both_strand_weak >= 0),
                (min_sd_both_strand_weak >= 0),
                (min_mad_both_strand_strong >= 0),
                (min_sd_both_strand_strong >= 0),
                (min_reads >= 0)]):
        logging.error('extended arg out of range, check helptext for ranges')
        sys.exit(EXIT_FAILURE)

    # set up filter params based on inputs TODO: move validation there, rather than above validation check
    al_params = ALF.Params(
        al_filter_threshold
    )
    ad_params = ADF.Params(
        edge_definition,
        edge_fraction,
        min_mad_one_strand,
        min_sd_one_strand,
        min_mad_both_strand_weak,
        min_sd_both_strand_weak,
        min_mad_both_strand_strong,
        min_sd_both_strand_strong,
        min_reads
    )
    dv_params = DVF.Params(
        duplication_window_size
    )

    try:
        vcf_in_handle = pysam.VariantFile(input_vcf_path)
    except Exception as er:
        logging.error(f'failed to open input VCF, reporting {er}')
        sys.exit(EXIT_FAILURE)
    vcf_names: list[str] = list(vcf_in_handle.header.samples)
    if len(set(vcf_names)) != len(vcf_names):
        logging.error('duplicate sample names in input VCF')
        sys.exit(EXIT_FAILURE)

    sm_to_aln_map: dict[str, pysam.AlignmentFile] = {}
    match format:
        case "s":
            mode = "r"
            logging.info("SAM format specified")
        case "b":
            mode = "rb"
            logging.info("BAM format specified")
        case "c":
            mode = "rc"
            logging.info("CRAM format specified")
        case _:
            raise ValueError(f'Unknown mode {format} found in format')
    for path in alignment_paths:
        try:
            alignment = pysam.AlignmentFile(path,
                                            mode,
                                            reference_filename=(cram_reference_path
                                                                if cram_reference_path
                                                                and format == "c"
                                                                else None))
        except Exception as er:
            logging.error(f'failed to read alignment file {path!r}, reporting {er}')
            sys.exit(EXIT_FAILURE)
        # grab the sample name from first SM field
        # in header field RG
        aln_sm = alignment.header.to_dict()['RG'][0]['SM']
        sm_to_aln_map[aln_sm] = alignment

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}
    if name_mapping:
        vcf_mapflag: list[str] = []
        alignment_mapflag: list[str] = []
        if len(name_mapping) <= len(alignment_paths) and all(m.count(':') == 1 for m in name_mapping) and not any("," in m for m in name_mapping):
            for pair in name_mapping:
                kv_split = pair.split(':')  # VCF:aln
                vcf_mapflag.append(kv_split[0])
                alignment_mapflag.append(kv_split[1])

            if len(vcf_mapflag) != len(set(vcf_mapflag)):
                logging.error('duplicate VCF sample names provided to --name-mapping flag')
                sys.exit(EXIT_FAILURE)

            if not set(vcf_mapflag) <= set(vcf_names):
                logging.error(f'VCF sample names provided to --name-mapping flag {vcf_mapflag!r} are not equal to or a subset of VCF samples from input VCF {vcf_names!r}')
                sys.exit(EXIT_FAILURE)

            if len(alignment_mapflag) != len(set(alignment_mapflag)):
                logging.error(msg='duplicate aligment sample names provided to name mapping flag')
                sys.exit(EXIT_FAILURE)

            if not sorted(alignment_mapflag) == sorted(sm_to_aln_map.keys()):  # check equal
                logging.error(msg='SM tags provided to name mapping flag {} are not equal to SM tag list from alignment files {} - flag recieved {}'.format(
                        alignment_mapflag,
                        set(sm_to_aln_map.keys()),
                        name_mapping))
                sys.exit(EXIT_FAILURE)

            vcf_sample_to_alignment_map = {vcf_mapflag[alignment_mapflag.index(k)]: v
                                            for k, v
                                            in sm_to_aln_map.items()}
        elif not any(":" in m for m in name_mapping) and len(alignment_paths) == len(name_mapping) == 1:
            if "," in name_mapping[0]:
                possible_sample_names = name_mapping[0].split(',')
            else:
                possible_sample_names = name_mapping  # list of len 1
            matches = [n for n in possible_sample_names if n in set(vcf_names)]
            if len(matches) > 1:
                logging.error(msg='More than one of the VCF sample names provided to name mapping flag match any sample names in input VCF {} - flag recieved {}'.format(
                        vcf_names,
                        name_mapping))
                sys.exit(EXIT_FAILURE)
            elif not matches:
                logging.error(msg='None of the VCF sample names provided to name mapping flag match any sample names in input VCF {} - flag recieved {}'.format(
                        vcf_names,
                        name_mapping))
                sys.exit(EXIT_FAILURE)
            else:
                logging.info('matched alignment to sample {} in VCF'.format(matches[0]))
                vcf_sample_to_alignment_map[matches[0]] = alignment  # pyright: ignore[reportPossiblyUnboundVariable] | since length of alignments == 1, can just reuse this variable

        else:
            logging.error(msg='name mapping misformatted, see helptext for expectations - flag recieved: {}'.format(name_mapping))
            sys.exit(EXIT_FAILURE)
    else:  # no name mapping flag
        if not sm_to_aln_map.keys() <= set(vcf_names):
            logging.error(
                msg='alignment SM tags {} are not equal to or a subset of VCF sample names {}'.format(
                    set(sm_to_aln_map.keys()),
                    set(vcf_names)
                )
            )
            sys.exit(EXIT_FAILURE)
        vcf_sample_to_alignment_map = sm_to_aln_map
    ## end name mapping handling

    if set(vcf_names) != vcf_sample_to_alignment_map.keys():
        logging.info(
            "alignments not provided for all VCF samples; {} will be ignored".format(
                set(vcf_names) - vcf_sample_to_alignment_map.keys()
            )
        )

    # init output  TODO: put these in filter modules as funcs
    out_head = vcf_in_handle.header
    out_head.add_line(
        f'##FILTER=<ID=ALF,Description="Median alignment score of reads reporting variant less than {al_filter_threshold}, using samples {vcf_sample_to_alignment_map.keys()!r}">'
    )
    out_head.add_line(
        f'##FILTER=<ID=ADF,Description="Variant arises from hairpin artefact, using samples {vcf_sample_to_alignment_map.keys()}">'
    )
    out_head.add_line(
        '##FILTER=<ID=DVF,Description="PLACEHOLDER">'
    )  # TODO: placeholder!
    out_head.add_line(
        '##INFO=<ID=ADF,Number=1,Type=String,Description="alt|code for each alt indicating hairpin filter decision code">'
    )
    out_head.add_line(
        '##INFO=<ID=ALF,Number=1,Type=String,Description="alt|code|score for each alt indicating AL filter conditions">'
    )
    out_head.add_line(
        '##INFO=<ID=DVF,Number=1,Type=String,Description="alt|code|score for each alt indicating DV filter conditions">'
    )  # TODO: placeholder!
    out_head.add_line(
        f'##hairpin2_version={__version__}'
    )
    out_head.add_line(
        f'##hairpin2_params=[{json.dumps({k: arg_d[k] for k in rec_args})}]'
    )

    try:
        vcf_out_handle = pysam.VariantFile(sys.stdout, 'w', header=out_head)  # TODO: check stdout works here
    except Exception as e:
        logging.error(msg='failed to open VCF output, reporting: {}'.format(e))
        sys.exit(EXIT_FAILURE)

    # test records
    prog_bar_counter = 0
    for record in vcf_in_handle.fetch():
        record_filters: dict[type[FilterResult[Any]], list[FilterResult[Any]]] = {ALF.Result: [], ADF.Result: [], DVF.Result: []}
        if record.alts is None:
            logging.warning('{0: <7}:{1: >12} ¦ no alts for this record, skipping'.format(
                record.chrom, record.pos))
        else:
            samples_w_mutants = [name
                                 for name
                                 in record.samples
                                 if record.samples[name]["GT"] != (0, 0)]
            if len(samples_w_mutants) == 0:
                logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit alts associated with this record, skipping'.format(
                    record.chrom, record.pos))
            else:
                region_reads_by_sample: dict[str, Iterable[pysam.AlignedSegment]] = {}
                for k, v in vcf_sample_to_alignment_map.items():
                    if k in samples_w_mutants:
                        read_iter, test_iter = tee(v.fetch(record.chrom,
                                                           record.start,
                                                           (record.start + 1)))
                        try:
                            _ = next(test_iter)
                        except StopIteration:
                            continue
                        else:
                            broad_qc_region_reads = [
                                read for read in read_iter
                                if not qc_read_broad(
                                    read,
                                    record.start,
                                    min_mapping_quality,
                                    min_clip_quality
                                )
                            ]
                            # doesn't check for overwrite
                            region_reads_by_sample[k] = broad_qc_region_reads

                # TODO: put mutation type detection under testing
                # the ability to handle complex mutations would be a potentially interesting future feature
                # for extending to more varied artifacts
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
                        logging.warning('could not infer mutation type, POS={} REF={} ALT={}, skipping variant'.format(
                            record.pos, record.ref, alt))
                        continue

                    alt_qc_region_reads: dict[str, list[pysam.AlignedSegment]] = {
                        key: [] for key in region_reads_by_sample
                    }
                    for mut_sample, read_iter in region_reads_by_sample.items():
                        alt_qc_region_reads[mut_sample] = [
                            read for read in read_iter
                            if not qc_read_alt_specific(
                                read,
                                record.start,
                                record.stop,
                                alt,
                                mut_type,
                                min_base_quality
                            )
                        ]

                    # instantiate filters to test the QC'd reads
                    ad = ADF.Filter(fixed_params=ad_params)
                    al = ALF.Filter(fixed_params=al_params)
                    dv = DVF.Filter(fixed_params=dv_params)

                    # run the tests
                    dup_qc_region_reads, dv_result = dv.test(alt, alt_qc_region_reads)
                    testing_reads = list(chain.from_iterable(dup_qc_region_reads.values()))
                    testing_reads, al_result = al.test(alt, testing_reads)
                    testing_reads, ad_result = ad.test(alt, record.start, testing_reads)
                    for res in (al_result, ad_result, dv_result):
                        record_filters[type(res)].append(res)

        if any(lst for lst in record_filters.values()):
            for ftype in record_filters:
                if any(fres.flag == True for fres in record_filters[ftype]):
                    record.filter.add(ftype.Name)
                record.info.update({ftype.Name: ','.join([fl.getinfo() for fl in record_filters[ftype]])})  # pyright: ignore[reportArgumentType]

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

    if output_json:
        try:
            with open(output_json, "w") as output_json:
                json.dump(
                    {
                        k: arg_d[k]
                        for k
                        in rec_args
                    },
                    output_json, indent="")
        except Exception as e:
            logging.error(msg='failed to write output JSON, reporting: {}'.format(e))
            sys.exit(EXIT_FAILURE)

    logging.info('hairpin complete')
    sys.exit(EXIT_SUCCESS)
