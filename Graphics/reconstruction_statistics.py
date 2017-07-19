#!/usr/bin/env python3

import argparse
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import scipy as sc
import numpy as np
from math import ceil, floor
from scipy.stats import hmean
from itertools import zip_longest
from collections import OrderedDict
from utils import parse_configuration
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import os
import re
# import matplotlib.gridspec as gridspec
# import matplotlib.ticker as ticker

__doc__ = """Script to automate the plots for the Mikado compare statistics"""


def calc_f1(rcl, prc):

    if rcl == 0 or prc == 0:
        return 0
    else:
        return 2*(rcl*prc)/(rcl+prc)


def split_comma(string):
    return string.split(",")


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def rcl(f1s, prc):
    denominator = (2.0 * prc - f1s)
    if denominator == 0:
        return 110
    else:
        return f1s * prc / denominator


def clamp(x):
    return max(0, min(x, 255))


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

def main():

    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--tight", default=False, action="store_true")
    parser.add_argument("--title", default="Mikado stats")
    parser.add_argument("--format", default=None, choices=["png", "pdf", "ps", "eps", "svg"])
    parser.add_argument("--levels", default=None, choices=["base", "exon", "intron",
                                                           "intron_chain", "transcript", "gene"],
                        nargs="+")
    parser.add_argument("--dpi", type=int, default=None)
    args = parser.parse_args()

    options = parse_configuration(args.configuration)
    options["out"] = os.path.splitext(args.out)[0]

    stats = OrderedDict()

    # figure = plt.figure(dpi=options["dpi"], figsize=(12, 6))

    # gs = gridspec.GridSpec()

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

    if args.levels is None:
        nrows = 2
        ncols = 3
        indices = [name_corr[_][1] for _ in name_corr]
    else:
        if len(set(args.levels)) <= 2:
            nrows = 1
            ncols = len(set(args.levels))
        elif len(set(args.levels)) == 3:
            ncols = 3
            nrows = 1
        elif len(set(args.levels)) == 4:
            nrows = 2
            ncols = 2
        else:
            nrows = 2
            ncols = 3
        assert nrows * ncols >= len(set(args.levels))
        __arr = []
        indices = []
        for key in name_corr:
            if key in args.levels:
                __arr.append(name_corr[key][0])
                indices.append(name_corr[key][1])
        name_ar = np.array(__arr)

    figure, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        dpi=options["dpi"],
        figsize=(10, 6))

    figure.suptitle(" ".join(["${}$".format(_) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")

    Xaxis = mpatches.FancyArrow(0.1, 0.19, 0.89, 0,
                                width=0.001,
                                length_includes_head=True,
                                transform=figure.transFigure, figure=figure,
                                color="k")
    Yaxis = mpatches.FancyArrow(0.1, 0.19, 0, 0.7,
                                width=0.001,
                                length_includes_head=True,
                                transform=figure.transFigure, figure=figure,
                                color="k")
    figure.lines.extend([Xaxis, Yaxis])

    figure.text(0.92, 0.21, "$Recall$", ha="center", fontsize=15)
    figure.text(0.07, 0.8, "$Precision$", va="center", fontsize=15, rotation="vertical")

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(options["colourmap"]["name"])
    # mapper = cm.ScalarMappable(colors, "PuOr")

    for xrow in range(nrows):
        for yrow in range(ncols):
            if nrows > 1:
                key = name_ar[xrow, yrow]
                plot = axes[xrow, yrow]
            else:
                key = name_ar[yrow]
                if ncols > 1:
                    plot = axes[yrow]
                else:
                    plot = axes
            if key is None:
                continue
            # plot.grid(True, linestyle='dotted')
            # plot.set(adjustable="box-forced", aspect="equal")
            plot.set_title("{} level".format(key), fontsize=10)

            # plot.set_xlabel("Precision", fontsize=10)
            # plot.set_ylabel("Recall", fontsize=10)
            stats[key] = dict()
            stats[key][b"plot"] = plot
            for division in options["divisions"]:
                stats[key][division.encode()] = []

    for method in options["methods"]:
        for aligner in options["divisions"]:
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

                stats[list(stats.keys())[index]][aligner.encode()].append((precision, recall, f1))



    divisions = sorted(options["divisions"].keys())

    handles = None
    best_marker = "o"

    for stat in stats.keys():

        plot = stats[stat][b"plot"]

        ys = [np.array([_[0] for _ in stats[stat][division.encode()]]) for division in divisions]
        xs = [np.array([_[1] for _ in stats[stat][division.encode()]]) for division in divisions]

        # print("Xs", xs)
        # print("Ys", ys)

        # plot.axis("scaled")
        # Select a suitable maximum
        index = np.argwhere(name_ar_orig == stat)
        index = index[0][0] * 3 + index[0][1]
        # print("Index", index, stat, name_ar)

        # Structure of X:
        # array[ method1[X assembler1, X assembler2, X assembler3 ...],
        #        method2[X assembler1, X assembler2, X assembler3 ...],
        #        ... ]

        suitable_x = []
        for method_num in range(len(xs)):
            for assembler_val in xs[method_num]:
                if assembler_val > -10:
                    suitable_x.append(assembler_val)

        suitable_y = []
        for method_num in range(len(ys)):
            for assembler_val in ys[method_num]:
                if assembler_val > -10:
                    suitable_y.append(assembler_val)

        suitable_x = np.array(suitable_x)
        suitable_y = np.array(suitable_y)

        # print(suitable_x.min(), index, suitable_x)

        if args.tight is True:
            margin = 1
        else:
            margin = 5
            
        x_minimum = max(0, floor(suitable_x.min()) - margin)
        y_minimum = max(0, floor(suitable_y.min()) - margin)

        x_maximum = min(100,
                        ceil(max(_.max() for _ in xs)) + margin)
        y_maximum = min(100,
                        ceil(max(_.max() for _ in ys)) + margin)

        plotf1curves(plot, fstepsize=ceil(min(x_maximum - x_minimum, y_maximum - y_minimum)/10))
        best_f1 = (-1, [])

        for enumerated, division in enumerate(divisions):
            for index, vals in enumerate(zip(xs[enumerated], ys[enumerated], options["methods"].keys())):
                x, y, label = vals
                f1 = calc_f1(x, y)
                if best_f1[0] < f1:
                    best_f1 = (f1, [(x, y)])
                elif best_f1[0] == f1:
                    best_f1[1].append((x, y))

                if options["colourmap"]["use"] is False:
                    colour = options["methods"][label]["colour"]
                    matched = re.match("\(([0-9]*), ([0-9]*), ([0-9]*)\)$", colour)
                    if matched:
                        colour = "#{0:02x}{1:02x}{2:02x}".format(clamp(int(matched.groups()[0])),
                                                                 clamp(int(matched.groups()[1])),
                                                                 clamp(int(matched.groups()[2])))
                elif options["methods"][label]["colour"] in ("black", "k"):
                    colour = "black"
                else:
                    colour = color_map(color_normalizer(options["methods"][label]["index"]))
                # The size has to change in a inversely related way compared to the number of plots
                # Top must be 100
                plot.scatter(x, y,
                             label=label,
                             # label="{0} ({1})".format(label, division),
                             c=colour, marker=options["divisions"][division]["marker"],
                             edgecolor="k", s=[150/max(nrows,ncols)], alpha=.8)

        circle_rad = 30
        for best in best_f1[1]:
            plot.plot(best[0], best[1], "o",
                      label="Best F1",
                      ms=circle_rad,
                      linestyle="-",
                      mec="k",
                      mfc="none")
        # for best in best_f1[1]:
        #     plot.scatter(best[0], best[1],
        #                  label="Best F1",
        #                  marker=best_marker, s=[20], c="k")

        if handles is None:
            handles, labels = plot.get_legend_handles_labels()

        # labels = list(divisions) + list(options["methods"].keys())

        __axes = plot.axes
        # print(stat, "({}, {})".format(x_minimum, x_maximum), "({}, {})".format(y_minimum, y_maximum))
        __axes.set_xlim(x_minimum, x_maximum)
        __axes.set_ylim(y_minimum, y_maximum)
        # __axes.set_aspect("equal")
        plot.tick_params(axis='both', which='major', labelsize=8)

    # Now create the labels
    # First the aligners

    # labels = []

    div_labels = []

    f1_line = mlines.Line2D([], [], color="gray", linestyle="--")
    div_labels.append((f1_line, "F1 contours"))

    for division in options["divisions"]:
        faux_line = mlines.Line2D([], [], color="white",
                                  marker=options["divisions"][division]["marker"],
                                  markersize=14,
                                  markerfacecolor="black")
        div_labels.append((faux_line, division))

    # best_marker_line = mlines.Line2D([], [], color="white",
    #                             marker=best_marker, markersize=6,
    #                             markerfacecolor="black", markeredgecolor="black")
    best_marker_line = mlines.Line2D([], [], color="white",
                                marker=best_marker, markersize=15,
                                markerfacecolor="none", markeredgecolor="black")
    div_labels.append((best_marker_line, "Best F1"))

    for method in options["methods"]:
        if options["colourmap"]["use"] is False:
            colour = options["methods"][method]["colour"]
            matched = re.match("\(([0-9]*), ([0-9]*), ([0-9]*)\)$", colour)
            if matched:
                colour = "#{0:02x}{1:02x}{2:02x}".format(clamp(int(matched.groups()[0])),
                                                         clamp(int(matched.groups()[1])),
                                                         clamp(int(matched.groups()[2])))
        elif options["methods"][label]["colour"] in ("black", "k"):
            colour = "black"
        else:
            colour = color_map(color_normalizer(options["methods"][method]["index"]))

        patch = mpatches.Patch(facecolor=colour, linewidth=1, edgecolor="k")
        div_labels.append((patch, method))

    npatches = len(options["divisions"]) + len(options["methods"]) + 2
    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="lower center",
                  scatterpoints=1,
                  ncol=ceil(npatches/4),
                  fontsize=10,
                  framealpha=0.5)
    # Necessary to pad the superior title
    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.1,  # Left
                           0.2,  # Bottom
                           0.85,  # Right
                           0.9])  # Top
    if options["out"] is None:
        plt.ion()
        plt.show(block=True)
    else:
        if args.format is not None:
            options["format"] = args.format
        if args.dpi is not None:
            options["dpi"] = args.dpi

        plt.savefig("{}.{}".format(options["out"], options["format"]),
                    format=options["format"],
                    dpi=options["dpi"],
                    transparent=options["opaque"])

if __name__ == "__main__":
    main()
