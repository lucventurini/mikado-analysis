#!/usr/bin/env python3

import argparse
import itertools
import matplotlib.lines
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
from collections import OrderedDict
from utils import parse_configuration
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import yaml
import re
import warnings
import sys
from operator import itemgetter
import os


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

    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # The species should be a configuration file listing
    # for each species the following:
    # - folder
    # - name
    parser.add_argument("--species", type=argparse.FileType("r"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--opaque", default=True, action="store_false")
    parser.add_argument("--dpi", default=1000, type=int)
    parser.add_argument("--format", default="svg", choices=["png", "pdf", "ps", "eps", "svg"])
    parser.add_argument("--level", default="transcript",
                        choices=["base", "exon", "intron", "intron_chain", "transcript", "gene"])
    parser.add_argument("-s", "--style", choices=plt.style.available, default="ggplot")
    parser.add_argument("--title", default="Mikado stats for transcript level")
    args = parser.parse_args()

    species = yaml.load(args.species)

    names = species.pop("names")

    keys = dict()
    for name in species:
        if isinstance(species[name], dict) and species[name].get("name", None) in names:
            keys[species[name]["name"]] = name

    ncols = len(names)

    figure, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        dpi=args.dpi,
        figsize=(8 / 2 * ncols, 4 / 2 * ncols))

    figure.suptitle(" ".join(["${}$".format(_) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")
    figure.text(0.5, 0.15, "${}\ \%$".format(args.level.title()), ha="center", fontsize=15)

    # name_ar = []
    # This is tricky .. this should be divided in 2 if there is more than one level

    # categories = species.pop("categories")
    # for category in categories:
    #     for name in names:
    #         key = [_ for _ in species.keys() if species[_]["name"] == name and species[_]["category"] == category]
    #         assert len(key) == 1
    #         key = key.pop()
    #         name_ar.append(key)
    # name_ar = np.array(list(grouper(name_ar, ceil(len(species) / 2), None)))

    # Dictionary to indicate which line should be taken given the level
    line_correspondence = {"base": 5,
                           "exon": 7,
                           "intron": 8,
                           "intron_chain": 9,
                           "transcript": 12,
                           "gene": 15}

    # best_marker = "o"

    for yrow, name in enumerate(names):
        plot = axes[yrow]
        plot.set_title("${}$".format(names[yrow]), fontsize=15)
        stats = OrderedDict()

        with open(species[keys[name]]["configuration"]) as configuration:
            options = parse_configuration(configuration, prefix=species[keys[name]]["folder"])

            # Usually STAR, TopHat
            method_names = sorted(list(itertools.product(options["methods"], options["divisions"])))

            for division in options["divisions"]:
                stats[division.encode()] = []

            max_point = float("-inf")
            for counter, (method, division) in enumerate(reversed(method_names), 1):
                    # print("Method:", method, "Aligner:", division)
                try:
                    orig, filtered = options["methods"][method][division]
                except TypeError:
                    warnings.warn("Something went wrong for {}, {}; continuing".format(
                        method, division))
                    stats[division.encode()].append((-10, -10, -10))
                    continue
                orig_lines = [line.rstrip() for line in open(orig)]
                filtered_lines = [line.rstrip() for line in open(filtered)]
                # for index, line_index in enumerate([5, 7, 8, 9, 11, 12, 14, 15]):

                plot.plot((0, 100), (counter, counter), c="lightgrey")

                for index, line_index in enumerate([line_correspondence[args.level]]):
                    precision = float(orig_lines[line_index].split(":")[1].split()[1])
                    recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                    try:
                        f1 = hmean(np.array([precision, recall]))
                    except TypeError as exc:
                        raise TypeError("\n".join([str(_) for _ in [(precision, type(precision)),
                                                                    (recall, type(recall)),
                                                                    exc]]))
                    # print(level, method, division, (precision, recall, f1))
                    # stats[division.encode()].append((precision, recall, f1))
                    # We can plot directly
                    max_point = max([max_point] + [_ + 7 for _ in (precision, recall, f1)])
                    plot.scatter(precision, counter,
                                 label="Precision",
                                 # label="{0} ({1})".format(label, division),
                                 c="orange", marker="o",
                                 edgecolor="k", s=[100.0], alpha=1)
                    plot.scatter(recall, counter,
                                 label="Recall",
                                 # label="{0} ({1})".format(label, division),
                                 c="lightblue", marker="o",
                                 edgecolor="k", s=[100.0], alpha=1)
                    plot.scatter(f1, counter,
                                 label="F1",
                                 # label="{0} ({1})".format(label, division),
                                 c="black", marker="o",
                                 edgecolor="k", s=[100.0], alpha=1)

            __axes = plot.axes
            max_point = min(100, max_point)
            print(name, max_point)
            __axes.set_xlim(0, max_point)
            __axes.set_ylim(0, len(method_names) + 1)
            __axes.spines["top"].set_visible(False)
            __axes.spines["right"].set_visible(False)
            # __axes.set_aspect("equal")
            if yrow == 0:
                __axes.set_yticks(range(1, len(method_names) + 1))
                __axes.set_yticklabels(["{} ({})".format(*_) for _ in reversed(method_names)],
                                       # fontdict= {'family': 'serif',
                                       #            'color': 'black',
                                       #            'weight': 'normal',
                                       #            'size': 200},
                                       fontsize=16
                                       )
            else:
                __axes.spines["left"].set_visible(False)
                __axes.set_yticks([])

            plot.tick_params(axis='x', which='major', labelsize=8)

    div_labels = []

    # f1_line = mlines.Line2D([], [], color="gray", linestyle="--")
    # div_labels.append((f1_line, "F1 contours"))

    prec_line = mlines.Line2D([], [],  marker="o", markersize=20, markerfacecolor="orange", color="white")
    rec_line = mlines.Line2D([], [],  marker="o", markersize=20, markerfacecolor="lightblue", color="white")
    f1_line = mlines.Line2D([], [],  marker="o", markersize=20, markerfacecolor="black", color="white")

    div_labels.append((prec_line, "Precision"))
    div_labels.append((rec_line, "Recall"))
    div_labels.append((f1_line, "F1"))

    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="upper center",
                  scatterpoints=1,
                  ncol=3,
                  fontsize=14,
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
