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
# @click.help_option('-h', '--help', help='show helptext')
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
    metavar= 'S:SM S:SM...',
    help="If sample names in VCF differ from SM tags in alignment files, provide a key here to map them. "
    "When multiple alignments are provided, accepts a space separated list of sample:SM pairs. "
    "When only a single alignment is provided, also accepts a comma separated string of one or more possible sample-of-interest "
    "names like TUMOR,TUMOUR"
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

# === Config override groups (hidden by default) ===
@optgroup.group("\nread validation config overrides", hidden=True)
@optgroup.option(
    '--min-clip-quality',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_clip_quality'].range),
    show_default=str(_PARAMS['min_clip_quality'].default),
    help='discard reads with mean base quality of aligned bases below this value, if they have soft-clipped bases'
)
@optgroup.option(
    '--min-mapping-quality',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_mapping_quality'].range),
    show_default=str(_PARAMS['min_mapping_quality'].default),
    help='discard reads with mapping quality below this value'
)
@optgroup.option(
    '--min-base-quality',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_base_quality'].range),
    show_default=str(_PARAMS['min_base_quality'].default),
    help='discard reads with base quality below this value at variant position'
)

@optgroup.group('\nDVF config overrides', hidden=True)
@optgroup.option(
    '--duplication-window-size',
    metavar='INT',
    type=click.IntRange(*_PARAMS['duplication_window_size'].range),
    show_default=str(_PARAMS['duplication_window_size'].default),
    help='inclusive maximum window size, in number of bases to use when detecting PCR duplicates. -1 will disable duplicate detection'
)  # TODO/BUG: is -1 still true?

@optgroup.group("\nALF config overrides", hidden=True)
@optgroup.option(
    '--al-filter-threshold',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['al_filter_threshold'].range),
    show_default=str(_PARAMS['al_filter_threshold'].default),
    help='threshold for median of read alignment score per base of all relevant reads, at and below which a variant is flagged as ALF'
)

@optgroup.group("\nADF config overrides", hidden=True)
@optgroup.option(
    '--edge-definition',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['edge_definition'].range),
    show_default=str(_PARAMS['edge_definition'].default),
    help='percentage of a read that is considered to be "the edge" for the purposes of assessing variant location distribution'
)
@optgroup.option(
    '--edge-fraction',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['edge_fraction'].range),
    show_default=str(_PARAMS['edge_fraction'].default),
    help='percentage of variants must occur within EDGE_FRACTION of read edges to mark ADF flag'
)
@optgroup.option(
    '--min-mad-one-strand',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_mad_one_strand'].range),
    show_default=str(_PARAMS['min_mad_one_strand'].default),
    help='min range of distances between variant position and read start for valid reads when only one strand has sufficient valid reads for testing'
)  # BUG: what on earth does this helptext mean?
@optgroup.option(
    '--min-sd-one-strand',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['min_sd_one_strand'].range),
    show_default=str(_PARAMS['min_sd_one_strand'].default),
    help='min stdev of variant position and read start for valid reads when only one strand has sufficient valid reads for testing'
) # BUG: exclusive or not?!
@optgroup.option(
    '--min-mad-both-strand-weak',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_mad_both_strand_weak'].range),
    show_default=str(_PARAMS['min_mad_both_strand_weak'].default),
    help='min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbsw is true'
)
@optgroup.option(
    '--min-sd-both-strand-weak',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['min_sd_both_strand_weak'].range),
    show_default=str(_PARAMS['min_sd_both_strand_weak'].default),
    help='min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbsw is true- default: 2, range: 0-, exclusive'
)
@optgroup.option(
    '--min-mad-both-strand-strong',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_mad_both_strand_strong'].range),
    show_default=str(_PARAMS['min_mad_both_strand_strong']),
    help='min range of distances between variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -sbss is true'
)
@optgroup.option(
    '--min-sd-both-strand-strong',
    metavar='FLOAT',
    type=click.FloatRange(*_PARAMS['min_sd_both_strand_strong'].range),
    show_default=str(_PARAMS['min_sd_both_strand_strong'].default),
    help='min stdev of variant position and read start for valid reads when both strands have sufficient valid reads for testing AND -mbss is true'
)
@optgroup.option(
    '--min-reads',
    metavar='INT',
    type=click.IntRange(*_PARAMS['min_reads'].range),
    show_default=_PARAMS['min_reads'].default,
    help='number of reads at and below which the hairpin filtering logic considers a strand to have insufficient reads for testing'
)
def hairpin2(
    vcf: str,
    alignments: list[str],
    config_path: str | None,
    output_config_path: str | None,
    name_mapping: str | None,
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

    # override from config and ensure all args are set
    # test args are sensible, exit if not
    # unfortunately necessary duplication with options since the config remains unverified
    # TODO: move to validation of FilterParams with pydantic, where there should be full validation anyway
    # TODO: forbid extra config params (pydantic?)
    for key in kwargs:
        if key not in _PARAMS:
            logging.error(f'unrecognised parameter {key!r} - check spelling?')
        if kwargs[key] is None:
            if key in configd.keys():
                kwargs[key] = configd[key]
            elif key in _PARAMS:
                if not quiet: logging.info(f'parameter {key!r} not found in config or present on command line, falling back to standard default {_PARAMS[key].default}')
                kwargs[key] = _PARAMS[key].default
        if key in _PARAMS:
            assert kwargs[key] is not None
            pmin, pmax = _PARAMS[key].range
            if (pmin is not None and pmin > kwargs[key]) or (pmax is not None and kwargs[key] > pmax):
                logging.error(f'value {kwargs[key]!r} for parameter {key!r} out of range {pmin}<=x<={pmax}')
                sys.exit(EXIT_FAILURE)

    # set up filter params based on inputs
    al_params = ALF.Params(
        kwargs['al_filter_threshold']
    )
    ad_params = ADF.Params(
        kwargs['edge_definition'],
        kwargs['edge_fraction'],
        kwargs['min_mad_one_strand'],
        kwargs['min_sd_one_strand'],
        kwargs['min_mad_both_strand_weak'],
        kwargs['min_sd_both_strand_weak'],
        kwargs['min_mad_both_strand_strong'],
        kwargs['min_sd_both_strand_strong'],
        kwargs['min_reads']
    )
    dv_params = DVF.Params(
        kwargs['duplication_window_size']
    )

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

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}
    if name_mapping:
        name_mapping: list[str] = name_mapping.split(' ')
        vcf_mapflag: list[str] = []
        alignment_mapflag: list[str] = []
        if len(name_mapping) <= len(alignments) and all(m.count(':') == 1 for m in name_mapping) and not any("," in m for m in name_mapping):
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
        elif not any(":" in m for m in name_mapping) and len(alignments) == len(name_mapping) == 1:
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
                if not quiet: logging.info('matched alignment to sample {} in VCF'.format(matches[0]))
                vcf_sample_to_alignment_map[matches[0]] = alignment  # pyright: ignore[reportPossiblyUnboundVariable] | since length of alignments == 1, can just reuse this variable

        else:
            breakpoint()
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
        if not quiet: logging.info(
            "alignments not provided for all VCF samples; {} will be ignored".format(
                set(vcf_names) - vcf_sample_to_alignment_map.keys()
            )
        )

    # init output  TODO: put these in filter modules as funcs
    out_head = vcf_in_handle.header
    out_head.add_line(
        f'##FILTER=<ID=ALF,Description="Median alignment score of reads reporting variant less than {kwargs['al_filter_threshold']}">'
    )
    out_head.add_line(
        f'##FILTER=<ID=ADF,Description="Variant shows anomalous distribution in supporting reads">'
    )
    out_head.add_line(
        '##FILTER=<ID=DVF,Description="Variant consistent with PCR stuttering on supporting reads">'
    )
    out_head.add_line(
        '##INFO=<ID=ADF,Number=1,Type=String,Description="alt|[True,False]|code indicating decision reason for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=ALF,Number=1,Type=String,Description="alt|[True,False]|code|score indicating decision reason and average AS of supporting reads for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=DVF,Number=1,Type=String,Description="alt|[True,False]|code indicating decision reason for each alt">'
    )  # TODO: placeholder!
    out_head.add_line(
        f'##hairpin2_version={__version__}'
    )
    out_head.add_line(
        f'##hairpin2_params=[{json.dumps({k: kwargs[k] for k in _PARAMS})}]'
    )
    out_head.add_line(
        f'##hairpin2_samples={vcf_sample_to_alignment_map.keys()}'
    )

    try:
        vcf_out_handle = pysam.VariantFile(sys.stdout, 'w', header=out_head)
    except Exception as e:
        logging.error(msg='failed to open VCF output, reporting: {}'.format(e))
        sys.exit(EXIT_FAILURE)

    # test records
    prog_bar_counter = 0
    for record in vcf_in_handle.fetch():
        record_filters: dict[type[FilterResult[Any]], list[FilterResult[Any]]] = {ALF.Result: [], ADF.Result: [], DVF.Result: []}
        if record.alts is None:
            if quiet < 1: logging.warning('{0: <7}:{1: >12} ¦ no alts for this record, skipping'.format(
                record.chrom, record.pos))
        else:
            samples_w_mutants = [name
                                 for name
                                 in record.samples
                                 if record.samples[name]["GT"] != (0, 0)]
            if len(samples_w_mutants) == 0:
                if quiet < 1: logging.warning('{0: <7}:{1: >12} ¦ no samples exhibit alts associated with this record, skipping'.format(
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
                                    kwargs['min_mapping_quality'],
                                    kwargs['min_clip_quality']
                                )
                            ]
                            # TODO: doesn't check for overwrite
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
                        if quiet < 1: logging.warning('could not infer mutation type, POS={} REF={} ALT={}, skipping variant'.format(
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
                                kwargs['min_base_quality']
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

    if output_config_path:
        try:
            with open(output_config_path, "w") as output_json:
                json.dump(
                    {
                        k: kwargs[k]
                        for k
                        in _PARAMS
                    },
                    output_json, indent="")
        except Exception as e:
            logging.error(msg='failed to write output JSON, reporting: {}'.format(e))
            sys.exit(EXIT_FAILURE)

    if not quiet: logging.info('hairpin complete')
    sys.exit(EXIT_SUCCESS)
