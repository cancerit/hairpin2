import base64
import json
import sys
import zlib

import click
import pysam

from hairpin2.sci_funcs import *


@click.group("hp2-utils")
def hp2_utils():
    pass


# NOTE: an API which makes this not hacky to do is a good api... (but later release)
# @hp2_utils.command('explain-var')
# @click.argument('INFO')
# def explain_variant(
#     info: str
# ):
#     kv = info.split("=")
#     flag_name, infopack = kv

#     if flag_name not in FlaggerNamespaces:
#         print(f"info field does not appear to be a hairpin2 flag (flag name was {flag_name})")
#         sys.exit(1)

#     infol = infopack.split('|')
#     alt = infol[0]
#     outcome = infol[1]
#     code = int(infol[2])
#     nreads = infol[3]
#     extras = None
#     if len(infol) > 4:
#         extras = infol[3:]

#     decomposed_code = None
#     if flag_name == "ADF":
#         decomposed_code = [member.name for member in AnomalousDistributionTest.ResultPack.Info if member.value & code]

#     print(flag_name)
#     print(outcome)
#     print(alt)
#     print(decomposed_code)
#     print(nreads)
#     print(extras)


@hp2_utils.command("get-params")
@click.argument("VCF")
def get_params(vcf: str):
    # param_line = None
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
                    json.dump(configd, fp=sys.stdout, indent="  ")
                except Exception as ex:
                    print(
                        f"Failed to decode (base85 -> zlib -> json) hairpin2 parameters from value in header key hairpin2_params. Got value {val} and reported error: {ex}"
                    )
                sys.exit(0)
        else:
            print("Could not find hairpin2_params key in VCF header")
            sys.exit(1)


if __name__ == "__main__":
    hp2_utils()
