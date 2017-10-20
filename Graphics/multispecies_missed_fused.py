#!/usr/bin/env python3

import argparse
import itertools
import matplotlib.lines as lines
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import matplotlib.gridspec as gridspec
import scipy as sc
import numpy as np
import matplotlib.ticker as ticker
from math import ceil, floor
from scipy.stats import hmean
from itertools import zip_longest
from collections import OrderedDict, defaultdict
from utils import parse_configuration
import operator
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import yaml
import re
import warnings
import sys
from operator import itemgetter
import os
import pandas
import csv


__doc__ = """Script to automate the plots for the Mikado compare statistics"""


def split_comma(string):
    return string.split(",")


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def calc_f1(rcl, prc):

    if rcl == 0 or prc == 0:
        return 0
    else:
        return 2*(rcl*prc)/(rcl+prc)


def rcl(f1s, prc):
    denominator = (2.0 * prc - f1s)
    if denominator == 0:
        return 100
    else:
        return f1s * prc / denominator


def plotf1curves(axis, fstepsize=1, stepsize=0.1):
    p = np.arange(0.0, 100.0, stepsize)[1:]
    for f in sc.arange(0.0, 100.0, fstepsize)[1:]:
        points = [(x, rcl(f, x)) for x in p if 0 < rcl(f, x) <= 100.0]
        if len(points) > 0:
            xs, ys = zip(*points)
            #print("F1 grid")
            #print(xs, ys)
            axis.plot(xs, ys, "--", color="gray", linewidth=0.5)  # , label=r"$f=%.1f$"%f) # exclude labels, for legend
            # bad hack:
            # gets the 10th last datapoint, from that goes a bit to the left, and a bit down
            #axis.annotate(r"$f=%.1f$" % f, xy=(xs[-10], ys[-10]), xytext=(xs[-10] - 0.05, ys[-10] - 0.035), size="small", color="gray")


def clamp(x):
    return max(0, min(x, 255))


def main():

    parser = argparse.ArgumentParser(__doc__)
    # parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--species", type=argparse.FileType("r"))
    parser.add_argument("--title", required=True)
    parser.add_argument("--log", action="store_true", default=False)
    parser.add_argument("--out", required=True)
    parser.add_argument("--type", required=True, choices=["missed", "fused", "detected", "reconstructed"])
    parser.add_argument("--format", default="svg", choices=["png", "pdf", "ps", "eps", "svg"])
    parser.add_argument("--dpi", default=1000, type=int)
    parser.add_argument("--opaque", default=True, action="store_false")
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    species = yaml.load(args.species)

    names = species.pop("names")

    name_keys = dict()
    for name in species:
        if isinstance(species[name], dict) and species[name].get("name", None) in names:
            name_keys[species[name]["name"]] = name

    ncols = len(names)

    figure, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        dpi=args.dpi,
        figsize=(10 / 2 * ncols, 13))

    figure.suptitle(" ".join(["${}$".format(_) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")
    figure.text(0.5, 0.15, "${}\ genes\ by\ method,\ in\ other\ methods$".format(args.type.title()), ha="center", fontsize=15)

    for yrow, name in enumerate(names):
        plot = axes[yrow]
        plot.set_title("${}$".format(names[yrow]), fontsize=15)
        data = defaultdict(dict)
        with open(species[name_keys[name]]["configuration"]) as configuration:
            options = parse_configuration(configuration, prefix=species[name_keys[name]]["folder"])
        # options = parse_configuration(args.configuration)

        for label in options["methods"]:
            for aligner in options["divisions"]:
                print(label, aligner)
                if options["methods"][label][aligner] is not None:
                    # data[(label, aligner)] = [set(), set(), set()]
                    orig_refmap = "{}.refmap".format(
                        re.sub(".stats$", "", options["methods"][label][aligner][0]))
                    orig_refmap = pandas.read_csv(orig_refmap, delimiter="\t")
                    orig_refmap = orig_refmap[["ref_gene", "best_ccode"]].drop_duplicates("ref_gene")
                    orig_refmap.columns = ["ref_gene", "orig_ccode"]
                    orig_refmap.orig_ccode.fillna("NA", inplace=True)
                    filtered_refmap = "{}.refmap".format(
                        re.sub(".stats$", "", options["methods"][label][aligner][1]))
                    filtered_refmap = pandas.read_csv(filtered_refmap, delimiter="\t")
                    filtered_refmap = filtered_refmap[["ref_gene", "best_ccode"]].drop_duplicates("ref_gene")
                    filtered_refmap.columns = ["ref_gene", "filtered_ccode"]
                    filtered_refmap.filtered_ccode.fillna("NA", inplace=True)
                    conc = pandas.merge(orig_refmap, filtered_refmap, how="outer")
                    for row in conc.itertuples():
                        if row.orig_ccode in ("=", "_"):
                            category = 0
                        elif row.orig_ccode[0] == "f":
                            category = 2
                        elif not row.filtered_ccode is np.nan and row.filtered_ccode in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                            category = 3
                        else:
                            category = 1

                        data[row.ref_gene][(label, aligner)] = category

                    # for num in range(3):
                    #     data[(label, aligner)][num] = len(data[(label, aligner)][num])
                # else:
                #     data[(label, aligner)] = [-1000] * 3

        # Now we have to create the numbers
        counts = OrderedDict()

        type_dict = {"missed": 3, "fused": 2, "detected": 1, "reconstructed": 0}

        to_search = type_dict[args.type]

        # This is probably VERY inefficient
        keys = list(reversed(sorted(itertools.product(options["methods"], options["divisions"]),
                                    key=operator.itemgetter(1))))

        max_val = 0
        colours = ["darkorange", "darkcyan", "purple", "lightcoral"]
        for rowno, key in enumerate(keys, 1):
            # Full, detected, fused, missed
            counts[key] = [0, 0, 0, 0]
            for gene in data:
                assert isinstance(data[gene], dict), (gene, type(data[gene]))
                assert key in data[gene], (gene, key, data[gene].keys())
                if data[gene][key] == to_search:
                    for pos in range(4):
                        # We have to segregate by division
                        if any(data[gene][_] == pos for _ in data[gene] if _ != key and _[1] == key[1]):
                            counts[key][pos] += 1
                            break

            left = 0
            for num, old, colour in zip(counts[key], [0] + counts[key][:-1], colours):
                print(rowno, num, colour, left)
                left += old
                plot.barh(rowno, num, color=colour, left=left, height=0.5, edgecolor="k")
            max_val = max(max_val, sum(counts[key]))
            # max_val = max(max_val, sum(counts[key]))

        plot.axes.set_xlim(0, max_val * 1.1)
        plot.axes.set_ylim(0, len(keys) + 0.5)

        # Now draw the lines ...

        plot.spines["top"].set_visible(False)
        plot.spines["right"].set_visible(False)

        if yrow == 0:
            plot.axes.set_yticks(range(1, len(keys) + 1))
            plot.axes.set_yticklabels(["{}".format(_[0]) for _ in keys],
                                      fontsize=14)
            x, y = np.array([[-2000, -0.5], [0, 0]])
            line = lines.Line2D(x, y, lw=1., color='k', alpha=1)
            line.set_clip_on(False)
            plot.add_line(line)
            for div_index, div in enumerate(reversed(sorted(options["divisions"])), 1):
                line_pos = len(options["methods"]) * div_index + 0.5
                text_pos = len(options["methods"]) * (div_index - 0.5) + 0.5
                x, y = np.array([[-2000, -0.5], [line_pos, line_pos]])
                line = lines.Line2D(x, y, lw=1., color='k', alpha=1)
                line.set_clip_on(False)
                plot.add_line(line)
                # ax.plot((-1000, max_val), (line_pos, line_pos), color="lightgrey")
                plot.text(-max_val * 0.75, text_pos, div, ha="center", fontsize=15, rotation="vertical")
        else:
            plot.axes.set_yticks([])

    # Now the legend ...

    div_labels = []
    full = mpatches.Patch(facecolor=colours[0], linewidth=1, edgecolor="k")
    detected = mpatches.Patch(facecolor=colours[1], linewidth=1, edgecolor="k")
    fused = mpatches.Patch(facecolor=colours[2], linewidth=1, edgecolor="k")
    missed = mpatches.Patch(facecolor=colours[3], linewidth=1, edgecolor="k")

    div_labels.append((full, "Fully reconstructed by at least another method"))
    div_labels.append((detected, "Detected by at least another method"))
    div_labels.append((fused, "Fused or missed by any other method"))
    div_labels.append((missed, "Missed by all other methods"))

    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="lower center",
                  scatterpoints=1,
                  ncol=2,
                  fontsize=20,
                  framealpha=0.5)

    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.1,  # Left
                           0.2,  # Bottom
                           0.85,  # Right
                           0.9])  # Top
    if args.out is None:
        plt.ion()
        plt.show(block=True)
    else:
        plt.savefig("{}.{}".format(args.out, args.format),
                    format=args.format,
                    dpi=args.dpi,
                    transparent=args.opaque)


if __name__ == "__main__":
    main()
