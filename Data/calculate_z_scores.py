import os
from collections import OrderedDict
from scipy.stats import hmean, zscore
from utils import parse_configuration
import argparse
import numpy as np
import sys


def main():

    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("-t", "--type", default="F1", choices=["F1", "Precision", "Recall"])
    parser.add_argument("-o", "--out", default=None)
    args = parser.parse_args()

    options = parse_configuration(args.configuration)
    if args.out is not None:
        options["out"] = os.path.splitext(args.out)[0]
    else:
        options["out"] = None

    stats = OrderedDict()

    name_ar = np.array([["Base", "Exon", "Intron"],
                        ["Intron chain", "Transcript", "Gene"]])
    name_ar_orig = name_ar.copy()

    name_corr = OrderedDict()
    name_corr["base"] = ("Base", 5)
    name_corr["exon"] = ("Exon", 7)
    name_corr["intron"] = ("Intron", 8)
    name_corr["intron_chain"] = ("Intron chain", 9)
    name_corr["transcript"] = ("Transcript", 12)
    name_corr["gene"] = ("Gene", 15)
    indices = [name_corr[_][1] for _ in name_corr]

    for key in name_corr:
        stats[key] = []

    key_map = []

    for method in options["methods"]:
        for aligner in options["divisions"]:
            key_map.append((method, aligner))
            if options["methods"][method][aligner] is not None:
                orig, filtered = options["methods"][method][aligner]
                orig_lines = [line.rstrip() for line in open(orig)]
                filtered_lines = [line.rstrip() for line in open(filtered)]
            else:
                orig_lines = None
                filtered_lines = None
                print("Aligner {} not found for method {}".format(aligner, method))
            # for index, line_index in enumerate([5, 7, 8, 9, 11, 12, 14, 15]):
            for index, line_index in enumerate(indices):
                if orig_lines is not None:
                    precision = float(orig_lines[line_index].split(":")[1].split()[1])
                    recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                    try:
                        f1 = hmean(np.array([precision, recall]))
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

                stats[list(stats.keys())[index]].append((precision, recall, f1))

    zscores = OrderedDict()
    for stat in stats:
        stats[stat] = np.array(stats[stat])
        zscores[stat] = np.array([zscore(stats[stat][:, 0]), zscore(stats[stat][:, 1]), zscore(stats[stat][:, 2])])

    if options["out"] is not None:
        out = open(options["out"], "wt")
    else:
        out = sys.stdout

    stat_index = ["Precision", "Recall", "F1"].index(args.type)

    print("Method", *stats.keys(), sep="\t", file=out)
    for index, key in enumerate(key_map):
        line = []
        line.append("{} ({})".format(*key))
        for stat in zscores:
            line.append(zscores[stat][stat_index][index])
        print(*line, sep="\t")

    if options["out"] is not None:
        out.close()

main()