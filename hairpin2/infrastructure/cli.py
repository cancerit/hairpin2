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

# pyright: reportImplicitStringConcatenation=false

import base64
import json
import logging
import sys
import tomllib
import zlib
from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast, override

import click
import pysam

from hairpin2.infrastructure.rascheduler import ConfigError
from hairpin2.main import hairpin2
from hairpin2.process_wrappers.default_exec import DEFAULT_EXEC_CONFIG
from hairpin2.VERSION import VERSION

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s ¦ %(levelname)-8s ¦ %(message)s", datefmt="%I:%M:%S"
)


EXIT_SUCCESS = 0
EXIT_FAILURE = 1


existing_file_path = click.Path(
    exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True, path_type=Path
)

writeable_file_path = click.Path(
    dir_okay=False,
    writable=True,
    resolve_path=True,
)


class JSONOrFile(click.ParamType):
    @override
    def convert(self, value: Any, param: Any, ctx: Any):
        # try to interpret as a path
        try:
            path = cast(Path, existing_file_path.convert(value, param, ctx))
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except click.BadParameter:
            pass  # not a valid file path, treat as raw JSON
        except Exception as e:
            self.fail(f"Failed to parse JSON file '{value}': {e}", param, ctx)

        # try parsing as raw JSON string
        try:
            return json.loads(value)
        except json.JSONDecodeError as e:
            self.fail(f"Invalid JSON input: {e}", param, ctx)


class ConfigFile(click.ParamType):
    name: str = "config-file"

    @override
    def convert(self, value: Any, param: Any, ctx: Any):
        path = cast(Path, existing_file_path.convert(value, param, ctx))
        data = path.read_bytes()
        ext = path.suffix.lower()
        try:
            if ext == ".toml":
                return tomllib.loads(data.decode("utf-8"))
            elif ext == ".json":
                return json.loads(data)
            else:
                self.fail("Expected .toml or .json", param, ctx)
        except Exception as ex:
            self.fail(f"Failed to parse {path.name}: {ex}", param, ctx)


def _coerce_value(path: str, raw: str) -> Any:
    s = raw.strip()

    try:
        return json.loads(s)
    except Exception:
        logging.warning(
            f'--set {path}: could not interpet value "{raw}" as JSON-format value, assuming string'
        )
        return s


def _set_dleaf(targetd: dict[str, Any], path: str, value: Any) -> None:
    """
    Set obj[path] = value where path is 'a.b.c' and only override existing leaf.
    """
    if not path or "." in (path[0], path[-1]):
        raise click.BadParameter(f"Invalid path for --set: {path!r}")

    dot_segments = path.split(".")
    if any(seg == "" for seg in dot_segments):
        raise click.BadParameter(f"Invalid path (empty segment) for --set: {path!r}")
    d_recursor: dict[str, Any] = targetd

    # walk all but last key
    for seg in dot_segments[:-1]:
        if not isinstance(d_recursor, dict):
            raise click.BadParameter(f"--set {path}: '{seg}' is not a mapping/dict in config")  # pyright: ignore[reportUnreachable]
        if seg not in d_recursor:
            raise click.BadParameter(f"--set {path}: missing key '{seg}' in config")
        d_recursor = d_recursor[seg]

    terminal_key = dot_segments[-1]
    if not isinstance(d_recursor, dict):
        raise click.BadParameter(f"--set {path}: parent of '{terminal_key}' is not a dict")  # pyright: ignore[reportUnreachable]
    if terminal_key not in d_recursor:
        raise click.BadParameter(f"--set {path}: missing key '{terminal_key}' in config")
    if isinstance(d_recursor[terminal_key], dict) or isinstance(d_recursor[terminal_key], list):
        raise click.BadParameter(
            f"--set {path}: final target is a mapping/dict or array/list; only leaf values may be overridden"
        )

    # override
    d_recursor[terminal_key] = value


def _apply_overrides_callback(
    ctx: click.Context,
    param: click.Option,  # pyright: ignore[reportUnusedParameter]
    pairs: tuple[tuple[str, str], ...],
):
    configd: dict[str, Any] | None = ctx.params.get("configd")
    if not isinstance(configd, dict):
        raise click.BadParameter("--set must occur after --config")

    for keypath, raw in pairs:
        _set_dleaf(configd, keypath, _coerce_value(keypath, raw))

    return None


def _resolve_configd_callback(ctx: Any, param: Any, values: Iterable[dict[str, Any]]):  # pyright: ignore[reportUnusedParameter]
    """
    Merge config files
    """
    # handle substitutions, merge
    # TODO: substitutions from vars, substitutions from cli flags
    merged: dict[str, Any] = {}

    for configd in values:
        for key, val in configd.items():
            if (
                key in merged
            ):  # TODO: allow splitting params key only across files as long as no sub keys overlap
                raise click.BadParameter(
                    f"-c/--config: top-level key {key} appears in more than one config. Top-level keys may not be split across config files"
                )
            merged[key] = val

    if not merged.get("exec"):
        merged["exec"] = DEFAULT_EXEC_CONFIG
        logging.info("execution flow configuration not provided; using defaults")

    return merged


@click.command(
    epilog="see documentation at https://github.com/cancerit/hairpin2 or at tool install location for further information",
    options_metavar="[--help] [OPTIONS]",
)
# TODO: add short help opt -h!
@click.version_option(VERSION, "-v", "--version", message="%(version)s")
@click.argument("vcf", type=existing_file_path, required=True)
@click.argument("alignments", nargs=-1, type=existing_file_path, required=True)
# Procedural options
@click.option(
    "-c",
    "--config",
    "configd",
    multiple=True,
    metavar="FILEPATH",
    type=ConfigFile(),
    help="path to config TOML/s or JSON/s from which processes and execution will be configured. May be provided multiple times; individual configs can be provided for each top level key (params, exec).",
    callback=_resolve_configd_callback,
)
@click.option(
    "-o",
    "--output-config",
    "output_config_path",
    metavar="FILEPATH",
    type=writeable_file_path,
    help="log run configuration back to a new JSON file.",
)
@click.option(
    "--set",
    nargs=2,
    multiple=True,
    type=(str, click.UNPROCESSED),
    expose_value=False,
    callback=_apply_overrides_callback,
    help="Override values in the supplied config. Uses dot paths for the key, and JSON format for the value, e.g., --set params.DVF.read_loss_threshold 0.6. Must be provided after --config. May be provided multiple times.",
)
@click.option(
    "-m",
    "--name-mapping",
    metavar="JSON_STRING | FILEPATH",
    type=JSONOrFile(),
    help="If sample names in VCF differ from SM tags in alignment files, provide a key here to map them. "
    "Accepts a path to a JSON file, or JSON-formatted string of key-value pairs where keys are sample names in the VCF "
    "and all values are either the SM tag or the filepath of the relevant alignment "
    '- e.g. \'{"sample0": "PDxxA", "sample1": "PDxxB"}\' or \'{"sample0": "A.bam", ...}\'. '
    "When only a single alignment is provided, also accepts a JSON-spec top-level array of possible sample of interest names "
    '- e.g. \'["TUMOR","TUMOUR"]\'. '
    "Note that when providing a JSON-formatted string at the command line you must single quote the string, and use only double quotes internally.",
)
@click.option(
    "-r",
    "--cram-reference",
    "cram_reference_path",
    metavar="FILEPATH",
    type=existing_file_path,
    help="path to FASTA format CRAM reference, overrides $REF_PATH and UR tags for CRAM alignments.",
)
@click.option(
    "-q",
    "--quiet",
    count=True,
    help="be quiet (-q to not log INFO level messages, -qq to additionally not log WARN).",
    default=False,
)
@click.option(
    "-p",
    "--progress",
    "progress_bar",
    is_flag=True,
    help="display progress bar on stderr during run.",
    default=False,
)
def hairpin2_cli(
    vcf: str,
    alignments: list[str],
    configd: dict[str, Any],
    output_config_path: str | None,
    name_mapping: dict[str, str] | list[str] | None,
    cram_reference_path: str | None,
    quiet: int,
    progress_bar: bool,
) -> None:
    """
    read-aware artefactual variant flagging algorithms. Flag variants in VCF using statistics calculated from supporting reads found in ALIGNMENTS and emit the flagged VCF to stdout.
    """
    try:
        vcf_in_handle = pysam.VariantFile(vcf)
    except Exception as er:
        logging.error(f"failed to open input VCF, reporting {er}")
        sys.exit(EXIT_FAILURE)
    vcf_names: list[str] = list(vcf_in_handle.header.samples)
    if len(set(vcf_names)) != len(vcf_names):
        logging.error("duplicate sample names in input VCF")
        sys.exit(EXIT_FAILURE)

    sm_to_aln_map: dict[str, pysam.AlignmentFile] = {}
    filename_to_aln_map: dict[str, pysam.AlignmentFile] = {}
    for path_str in alignments:
        path = Path(path_str)
        try:
            match path.suffix[1]:
                case "s" | "S":
                    mode = "r"
                case "b" | "B":
                    mode = "rb"
                case "c" | "C":
                    mode = "rc"
                case _:
                    raise ValueError
        except (IndexError, ValueError):
            logging.error(
                f"Could not infer alignment format from suffix {path.suffix!r} of file {path_str!r}"
            )
            sys.exit(EXIT_FAILURE)
        if cram_reference_path and mode != "rc":
            logging.error(
                f"CRAM reference provided, but alignment at {path_str!r} not inferred as cram from suffix {path.suffix!r}"
            )
            sys.exit(EXIT_FAILURE)
        try:
            alignment = pysam.AlignmentFile(
                str(path),
                mode,
                reference_filename=(cram_reference_path if cram_reference_path else None),
            )
        except Exception as er:
            logging.error(f"failed to read alignment file {path!r}, reporting {er}")
            sys.exit(EXIT_FAILURE)
        # grab the sample name from the first SM field
        # in header field RG
        aln_sm = alignment.header.to_dict()["RG"][0]["SM"]  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        sm_to_aln_map[aln_sm] = alignment
        filename_to_aln_map[path.name] = alignment

    vcf_sample_to_alignment_map: dict[str, pysam.AlignmentFile] = {}

    match name_mapping:
        case list():
            matches = [name for name in name_mapping if name in set(vcf_names)]
            if len(matches) > 1:
                logging.error(
                    msg="More than one of the VCF sample names provided to name mapping {name_mapping} match any sample names in input VCF {vcf_names}!"
                )
                sys.exit(EXIT_FAILURE)
            elif not matches:
                logging.error(
                    msg=f"None of VCF sample names provided to name mapping {name_mapping} match any sample name in input VCF {vcf_names}!"
                )
                sys.exit(EXIT_FAILURE)
            else:
                if not quiet:
                    logging.info("matched alignment to sample {} in VCF".format(matches[0]))
                vcf_sample_to_alignment_map[matches[0]] = alignment  # pyright: ignore[reportPossiblyUnboundVariable] | since length of alignments == 1, can just reuse this variable
        case dict():
            vcf_map_keys: list[str] = []
            alignment_map_values: list[str] = []
            for key, val in name_mapping.items():
                vcf_map_keys.append(key)
                alignment_map_values.append(val)

            if len(vcf_map_keys) != len(set(vcf_map_keys)):
                logging.error("duplicate keys (VCF sample names) provided to name mapping!")
                sys.exit(EXIT_FAILURE)

            if not set(vcf_map_keys) <= set(vcf_names):
                logging.error(
                    f"keys (VCF sample names) provided to name mapping {vcf_map_keys!r} are not equal to or a subset of VCF samples from input VCF {vcf_names!r}!"
                )
                sys.exit(EXIT_FAILURE)

            if len(alignment_map_values) != len(set(alignment_map_values)):
                logging.error(
                    msg="duplicate values (alignment SM tags/filenames) provided to name mapping!"
                )
                sys.exit(EXIT_FAILURE)

            match_sm = sorted(alignment_map_values) == sorted(sm_to_aln_map.keys())
            match_fn = sorted(alignment_map_values) == sorted(filename_to_aln_map.keys())
            if match_sm:
                vcf_sample_to_alignment_map = {
                    vcf_map_keys[alignment_map_values.index(k)]: v for k, v in sm_to_aln_map.items()
                }
            elif match_fn:
                vcf_sample_to_alignment_map = {
                    vcf_map_keys[alignment_map_values.index(k)]: v
                    for k, v in filename_to_aln_map.items()
                }
            else:
                logging.error(
                    msg=f"values provided to name mapping {alignment_map_values!r} are not equal to the either the SM tags or the filenames of the alignments provided"
                )
                sys.exit(EXIT_FAILURE)

            if match_sm and match_fn:
                logging.warning(
                    f"intention ambiguous - values provided to name mapping {alignment_map_values!r} are equal to both the SM tags and the filenames of alignments. Mapping against SM tags, not filenames"
                )
        case None:
            if not sm_to_aln_map.keys() <= set(vcf_names):
                logging.error(
                    msg="alignment SM tags {} are not equal to or a subset of VCF sample names {} - use a name mapping".format(
                        set(sm_to_aln_map.keys()), set(vcf_names)
                    )
                )
                sys.exit(EXIT_FAILURE)
            vcf_sample_to_alignment_map = sm_to_aln_map
        case _:
            logging.error(
                f"name mapping misformatted. Expected deserialised JSON, received {type(name_mapping)}, containing {name_mapping}"
            )  # pyright: ignore[reportUnreachable]
            sys.exit(EXIT_FAILURE)

    if set(vcf_names) != vcf_sample_to_alignment_map.keys():
        if not quiet:
            logging.info(
                "alignments not provided for all VCF samples; {} will be ignored".format(
                    set(vcf_names) - vcf_sample_to_alignment_map.keys()
                )
            )

    # this should move
    # init output  TODO: put these in filter modules as funcs
    # TODO: placeholder variable accession, super fragile I know
    out_head = vcf_in_handle.header
    out_head.add_line(
        f'##FILTER=<ID=ALF,Description="Median alignment score of reads reporting variant less than {configd["params"]["ALF"]["avg_AS_threshold"]}">'  # TODO: move to result or flagger classes
    )
    out_head.add_line(
        '##FILTER=<ID=ADF,Description="Variant shows anomalous distribution in supporting reads">'
    )
    out_head.add_line(
        f'##FILTER=<ID=DVF,Description="More than {configd["params"]["DVF"]["read_loss_threshold"]} of reads supporting this variant are considered PCR stutter duplicates">'
    )
    out_head.add_line(
        f'##FILTER=<ID=LQF,Description="More than {configd["params"]["LQF"]["read_loss_threshold"]} of reads supporting this variant are considered low quality">'
    )
    out_head.add_line(
        '##INFO=<ID=ADF,Number=A,Type=String,Description="alt|outcome|flag|nreads|strand indicating decision reason for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=ALF,Number=A,Type=String,Description="alt|outcome|flag|nreads|score indicating decision reason and average AS of supporting reads for each alt">'
    )
    out_head.add_line(
        '##INFO=<ID=DVF,Number=A,Type=String,Description="alt|outcome|flag|nreads|loss indicating decision reason for each alt and ratio of supporting reads suspected to be duplicates">'
    )
    out_head.add_line(
        '##INFO=<ID=LQF,Number=A,Type=String,Description="alt|outcome|flag|nreads|loss indicating decision reason for each alt and ratio of supporting reads suspected to be low quality">'
    )
    # store complete description of hp2 run in header
    out_head.add_line(f"##hairpin2_version={VERSION}")
    out_head.add_line(f"##hairpin2_params_info=The hairpin2_params key contains a reduced representation of the parameters used to run hairpin2 on this VCF. It is a JSON compressed with zlib and encoded with base85. You can retreive the JSON with hp2_utils get-params.")
    # squeeze params into a tiny representation so we don't clog the output VCF
    confpack = zlib.compress(
        json.dumps(
            {
                "params": {
                    k: v for k, v in configd["params"].items() if configd["exec"][k]["enable"]
                },
                "exec": configd["exec"],
            },
            separators=(",", ":"),  # minify json
        ).encode("ascii"),
        level=9,
    )
    enc_confpack = base64.b85encode(  # encode bytes to ascii so we can store as valid string in VCF
        confpack  # NOTE: base91 even shorter, but not in stdlib
    ).decode("ascii")
    out_head.add_line(f"##hairpin2_params={enc_confpack}")
    out_head.add_line(f"##hairpin2_samples={set(vcf_sample_to_alignment_map.keys())}")

    # test records
    prog_counter = 0
    try:
        for flagged_record in hairpin2(
            vcf_in_handle.fetch(), vcf_sample_to_alignment_map, configd, quiet
        ):
            if prog_counter == 0:  # a little awkward, but error before dumping out the vcf header
                try:
                    vcf_out_handle = pysam.VariantFile(sys.stdout, "w", header=out_head)
                except Exception as e:
                    logging.error(msg="failed to open VCF output, reporting: {}".format(e))
                    sys.exit(EXIT_FAILURE)
            try:
                _ = vcf_out_handle.write(flagged_record)  # pyright: ignore[reportPossiblyUnboundVariable]
            except Exception as e:
                logging.error(msg="failed to write to vcf, reporting: {}".format(e))
                sys.exit(EXIT_FAILURE)

            if progress_bar:
                if not prog_counter % 100:
                    print(
                        f"\rchr: {flagged_record.chrom}, pos: {flagged_record.pos}",
                        end="",
                        flush=True,
                        file=sys.stderr,
                    )
            prog_counter += 1
        else:  # after exhausting records, print a line break
            if progress_bar:
                print(file=sys.stderr)
    except ConfigError as er:
        logging.error(f"failed to start hairpin2 executor from config/s, reporting: {er}")
        sys.exit(EXIT_FAILURE)

    if output_config_path:
        try:
            with open(output_config_path, "w") as output_json:
                json.dump(configd, output_json, sort_keys=False, indent="  ")
        except Exception as e:
            logging.error(msg="failed to write output JSON, reporting: {}".format(e))
            sys.exit(EXIT_FAILURE)

    if not quiet:
        logging.info("hairpin complete")
    sys.exit(EXIT_SUCCESS)


if __name__ == "__main__":
    hairpin2_cli()
