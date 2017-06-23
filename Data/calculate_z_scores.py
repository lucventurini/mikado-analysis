#!/usr/bin/env python3

import os
import operator
from collections import OrderedDict, defaultdict
from scipy.stats import hmean
from utils import parse_configuration
import argparse
import numpy as np
import scipy.stats as mstats
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

    stats = dict()

    for key in name_corr:
        stats[key] = OrderedDict()

    key_map = []

    for method in sorted(options["methods"].keys()):
        for aligner in sorted(options["divisions"].keys()):
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

                stats[list(stats.keys())[index]][(method, aligner)] = (precision, recall, f1)

    zscores = {"precision": defaultdict(OrderedDict),
               "recall": defaultdict(OrderedDict),
               "f1": defaultdict(OrderedDict)}

    for stat in stats:
        for index, feat in enumerate(["precision", "recall", "f1"]):
            vals = np.array([_[index] for _ in stats[stat].values()])
            for key, score in zip(stats[stat].keys(), mstats.zscore(vals)):
                zscores[feat][key][stat] = score

    # Now we have calculated the Zscore for each of the methods for each of the stats
    # Time to sum them up

    zscore_sum = defaultdict(OrderedDict)
    ranks = defaultdict(OrderedDict)
    for feat in ["precision", "recall", "f1"]:
        for key in zscores[feat]:
            zscore_sum[feat][key] = sum(zscores[feat][key].values())
        feat_ranks = sorted(zscore_sum[feat].items(), key=operator.itemgetter(1), reverse=True)
        feat_ranks = list(enumerate(feat_ranks, 1))
        d = dict()
        for rank, (key, zsum) in feat_ranks:
            d[key] = (zsum, rank)
        for key in stats["base"].keys():
            ranks[feat][key] = d[key][:]

    if options["out"] is not None:
        out = open(options["out"], "wt")
    else:
        out = sys.stdout

    print("Method", "Precision", "", "Recall", "", "F1", sep="\t", file=out)
    print("", *["zscore", "rank"] * 3, sep="\t", file=out)
    for key in ranks["precision"].keys():
        row = ["{} ({})".format(*key)]
        for feat in ["precision", "recall", "f1"]:
            row.extend(ranks[feat][key])
        print(*row, sep="\t", file=out)

    if options["out"] is not None:
        out.close()

main()