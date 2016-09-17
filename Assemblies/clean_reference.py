#!/usr/bin/env python3

from Mikado.parsers.GTF import GTF
import argparse
import yaml
import sys


__doc__ = """
Script to filter the human models from the GTF according to the YAML table of
acceptable biotypes for genes and transcripts.
"""


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("gtf")
    parser.add_argument("yaml")
    parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout, nargs="?")
    args = parser.parse_args()

    with open(args.yaml) as _:
        acceptable = yaml.load(_)
    assert isinstance(acceptable, dict)

    with open(args.gtf) as gtf:
        for row in GTF(args.gtf):
            if acceptable[row.attributes["gene_biotype"]][row.source] is True:
                print(row, file=args.out)
            else:
                continue
                continue
    return

main()
