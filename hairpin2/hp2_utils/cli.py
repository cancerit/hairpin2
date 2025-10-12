import base64
import json
import sys
import zlib

import click
import pysam

from hairpin2.const import FlaggerNamespaces
from hairpin2.sci_funcs import *


@click.group("hp2-utils")
def hp2_utils():
    pass


# thrown together at the last minute
# NOTE: an API which makes this not hacky to do is a good api... (but later release)
@hp2_utils.command("explain-var", short_help="explain a flagging decision")
@click.argument("INFO")
def explain_variant(info: str):
    """
    Given a hairpin2 info string from a VCF like "FLAG=ALT|OUTCOME|CONDITION|NREADS|...", explain the flagging decision.

    Returns result on stdout.
    """
    kv = info.split("=")
    flag_name, infopack = kv

    if flag_name not in FlaggerNamespaces:
        print(f"info field does not appear to be a hairpin2 flag (flag name was {flag_name})")
        sys.exit(1)

    infol = infopack.split("|")
    alt = infol[0]
    outcome = infol[1]
    code = int(infol[2], 16)
    nreads = infol[3]
    extras = None
    if len(infol) > 4:
        extras = infol[4:]
    assert extras

    decomposed_code = []
    match flag_name:
        case FlaggerNamespaces.ANOMALOUS_DISTRIBUTION:
            decomposed_code = [
                member.name
                for member in AnomalousDistributionTest.ResultPack.Info
                if member.value & code
            ]
            outcome_str = f"variant outcome was {outcome} via conditions {decomposed_code} on strand {extras[0]}"
        case FlaggerNamespaces.DUPLICATION:
            decomposed_code = [
                member.name for member in ProportionBasedTest.ResultPack.Info if member.value & code
            ]
            outcome_str = f"variant outcome was {outcome} via conditions {decomposed_code} with read loss of {extras[0]}"
        case FlaggerNamespaces.LOW_QUAL:
            decomposed_code = [
                member.name for member in ProportionBasedTest.ResultPack.Info if member.value & code
            ]
            outcome_str = f"variant outcome was {outcome} via conditions {decomposed_code} with read loss of {extras[0]}"
        case FlaggerNamespaces.POOR_ALIGNMENT_SCORE:
            decomposed_code = [
                member.name for member in AlignmentScoreTest.ResultPack.Info if member.value & code
            ]
            outcome_str = f"variant outcome was {outcome} via conditions {decomposed_code} with average alignment score of {extras[0]}"
        case _:
            print(f"info field does not appear to be a hairpin2 flag (flag name was {flag_name})")
            sys.exit(1)

    print(f"FLAG: {flag_name}")
    print(f"ALT: {alt}")
    print(outcome_str)
    print(f"reads examined: {nreads}")


@hp2_utils.command("get-params", short_help="get run parameters from a VCF")
@click.argument("vcf", metavar="VCF-PATH")
@click.option("--exec", is_flag=True)
def get_params(vcf: str, exec: bool = False):
    """
    Find and decode a hairpin2 parameter string from a VCF back into a JSON.

    Returns decoded JSON on stdout.
    """
    with pysam.VariantFile(vcf) as vf:
        for hrec in vf.header.records:
            if hrec.key == "hairpin2_params":
                val = hrec.value
                if val is None:
                    print(
                        "hairpin2_params key in VCF header does not appear to have a value associated"
                    )
                    sys.exit(1)
                try:
                    comp_bytes = base64.b85decode(val)
                    configd = json.loads(zlib.decompress(comp_bytes))
                    if not exec:
                        configd.pop('exec', None)
                    json.dump(configd, fp=sys.stdout, indent="  ")
                except Exception as ex:
                    print(
                        f"Failed to decode (base85 -> zlib -> json) hairpin2 parameters from value in header key hairpin2_params. Got value {val} and reported error: {ex}"
                    )
                return
        else:
            print("Could not find hairpin2_params key in VCF header")
            sys.exit(1)


if __name__ == "__main__":
    hp2_utils()
