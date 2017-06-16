#!/usr/bin/env python3

from utils import parse_configuration
import argparse
from scipy.stats import hmean
from math import ceil, floor
import sys
import numpy as np
from collections import OrderedDict


__doc__ = """This script allows to parse multiple STATs files and
report a collapsed table of recall, precision and F1."""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-c", "--conf", required=True)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    line_correspondence = {"base": 5,
                           "exon": 7,
                           "intron": 8,
                           "intron_chain": 9,
                           "transcript": 12,
                           "gene": 15}

    with open(args.conf) as configuration:
        options = parse_configuration(configuration)

    stats = OrderedDict()

    name_ar = np.array([["Base", "Exon", "Intron"],
                        ["Intron chain", "Transcript", "Gene"]])
    # name_ar_orig = name_ar.copy()

    name_corr = OrderedDict()
    name_corr["base"] = ("Base", 5)
    name_corr["exon"] = ("Exon", 7)
    name_corr["intron"] = ("Intron", 8)
    name_corr["intron_chain"] = ("Intron chain", 9)
    name_corr["transcript"] = ("Transcript", 12)
    name_corr["gene"] = ("Gene", 15)

    for key in name_corr:
        stats[key] = OrderedDict()
        for division in options["divisions"]:
            stats[key][division.encode()] = dict()

    for method in options["methods"]:
        for division in options["divisions"]:
            if options["methods"][method][division] is not None:
                orig, filtered = options["methods"][method][division]
                orig_lines = [line.rstrip() for line in open(orig)]
                filtered_lines = [line.rstrip() for line in open(filtered)]
            else:
                orig_lines = None
                filtered_lines = None
                print("Aligner {} not found for method {}".format(division, method))
            # for index, line_index in enumerate([5, 7, 8, 9, 11, 12, 14, 15]):
            for index, line_index in enumerate([name_corr[_][1] for _ in name_corr]):
                if orig_lines is not None:
                    precision = orig_lines[line_index]
                    precision = precision.split(":")

                    precision = float(orig_lines[line_index].split(":")[1].split()[1])
                    recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                    try:
                        f1 = round(hmean(np.array([precision, recall])), 2)
                    except TypeError as exc:
                        raise TypeError("\n".join([str(_) for _ in [(precision, type(precision)),
                                                                    (recall, type(recall)),
                                                                    exc]]))
                else:
                    precision = -10
                    recall = -10
                    f1 = -10

                # In order:
                # Name of the statistic:Base, Exon, etc
                # Name of the method
                stats[list(stats.keys())[index]][division.encode()][method.encode()] = (precision, recall, f1)

    first_row = ["Level"]
    for name in name_corr:
        first_row.extend(["", name_corr[name][0], ""])
    print(*first_row, sep="\t", file=args.out)

    second_row = [""] + ["Precision", "Recall", "F1"] * (int((len(first_row) -1 )/ 3))
    print(*second_row, sep="\t", file=args.out)

    for method in options["methods"]:
        for division in options["divisions"]:
            row = ["{} ({})".format(method, division)]
            for name in name_corr:
                row.extend(stats[name][division.encode()][method.encode()])
            print(*row, file=args.out, sep="\t")

    return

main()
