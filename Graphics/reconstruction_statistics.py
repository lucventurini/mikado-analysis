#!/usr/bin/env python3

import argparse
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import matplotlib.ticker as ticker
from math import ceil, floor
from itertools import zip_longest
from collections import OrderedDict
from utils import parse_configuration

__doc__ = """Script to automate the plots for the Mikado compare statistics"""


def split_comma(string):
    return string.split(",")


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def main():

    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--out", required=True)
    parser.add_argument("--title", default="Mikado stats")
    args = parser.parse_args()

    options = parse_configuration(args)
    options["out"] = args.out

    stats = OrderedDict()
    figure, axes = plt.subplots(nrows=2,
                                ncols=4,
                                dpi=options["dpi"], figsize=(12, 6))
    figure.suptitle(" ".join(["${}$".format(_) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")

    Xaxis = matplotlib.patches.FancyArrow(0.1, 0.19, 0.89, 0,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")
    Yaxis = matplotlib.patches.FancyArrow(0.1, 0.19, 0, 0.7,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")
    figure.lines.extend([Xaxis, Yaxis])

    figure.text(0.92, 0.21, "$Precision$", ha="center", fontsize=15)
    figure.text(0.07, 0.8, "$Recall$", va="center", fontsize=15, rotation="vertical")

    name_ar = np.array([["Base", "Exon", "Intron", "Intron chain"],
                        ["Transcript (95% nF1)", "Transcript (80% nF1)",
                         "Gene (95% nF1)", "Gene (80% nF1)"]])

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(options["colourmap"]["name"])
    # mapper = cm.ScalarMappable(colors, "PuOr")

    for xrow in (0, 1):
        for yrow in (0, 1, 2, 3):
            key = name_ar[xrow, yrow]
            plot = axes[xrow, yrow]
            plot.grid(True)
            plot.set_title("{} level".format(key), fontsize=10)

            # plot.set_xlabel("Precision", fontsize=10)
            # plot.set_ylabel("Recall", fontsize=10)
            stats[key] = dict()
            stats[key][b"plot"] = plot
            stats[key][b"STAR"] = []
            stats[key][b"TopHat"] = []

    for method in options["methods"]:
        for aligner in ("STAR", "TopHat"):
            orig, filtered = options["methods"][method][aligner]
            orig_lines = [line.rstrip() for line in open(orig)]
            filtered_lines = [line.rstrip() for line in open(filtered)]
            for index, line_index in enumerate([5, 7, 8, 9, 11, 12, 14, 15]):
                precision = float(orig_lines[line_index].split(":")[1].split()[1])
                recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                stats[
                    # Name of the statistic:Base, Exon, etc
                    list(stats.keys())[index]][aligner.encode()].append((precision, recall))

    handles, labels = None, None

    for stat in stats.keys():
        star_values = stats[stat][b"STAR"]
        th_values = stats[stat][b"TopHat"]
        plot = stats[stat][b"plot"]

        Xtop = np.array([_[0] for _ in th_values])
        Ytop = np.array([_[1] for _ in th_values])

        Xstar = np.array([_[0] for _ in star_values])
        Ystar = np.array([_[1] for _ in star_values])

        # plot.axis("scaled")
        # Select a suitable maximum

        x_minimum = max(0, 10 * (floor(min(Xtop.min(), Xstar.min()) / 10.0) - 0.5))
        y_minimum = max(0, 10 * (floor(min(Ytop.min(), Ystar.min()) / 10.0) - 0.5))

        x_maximum = min(100, 10 * (ceil(max(Xtop.max(), Xstar.max()) / 10.0) + 0.5))
        y_maximum = min(100, 10 * (ceil(max(Ytop.max(), Ystar.max()) / 10.0) + 0.5))

        x_maximum, y_maximum = [max(x_maximum, y_maximum)] * 2
        x_minimum, y_minimum = [min(x_minimum, y_minimum)] * 2

        # minimum = max(0, 10 * (floor(
        #         min(Xtop.min(), Xstar.min(), Ystar.min(), Ytop.min()
        #         ) / 10.0) - 0.5))
        # maximum = min(100, 10 *(ceil(
        #     max(Xtop.max(), Xstar.max(), Ystar.max(), Ytop.max()) / 10.0) + 0.5))
        # maximum = 100

        plot.set_xlim(x_minimum, x_maximum)
        plot.set_ylim(y_minimum, y_maximum)
        plot.tick_params(axis='both', which='major', labelsize=8)
        __axes = plot.axes
        __axes.xaxis.set_major_locator(ticker.MultipleLocator(
            ceil((x_maximum -x_minimum)/ 20) * 5))
        __axes.xaxis.set_minor_locator(ticker.MultipleLocator(ceil((x_maximum -x_minimum)/ 100) * 5))
        __axes.yaxis.set_major_locator(ticker.MultipleLocator(ceil((y_maximum - y_minimum)/ 20) * 5))
        __axes.yaxis.set_minor_locator(ticker.MultipleLocator(ceil((y_maximum - y_minimum)/ 100) * 5))

        plot.plot(0, 1, ls="--", c=".3")

        for index, vals in enumerate(zip(Xtop, Ytop, options["methods"].keys())):
            x, y, label = vals
            if options["colourmap"]["use"] is False:
                colour = options["methods"][label]["colour"]
            else:
                colour = color_map(color_normalizer(options["methods"][label]["index"]))

            plot.scatter(x, y, label="{0} (TopHat)".format(label), c=colour,
                         marker="^", edgecolor="k",
                         s=[50.0], alpha=.8)

        for index, vals in enumerate(zip(Xstar, Ystar, options["methods"].keys())):
            x, y, label = vals
            if options["colourmap"]["use"] is False:
                colour = options["methods"][label]["colour"]
            else:
                colour = color_map(color_normalizer(options["methods"][label]["index"]))
            plot.scatter(x, y, label="{0} (STAR)".format(label), c=colour,
                         marker="o",
                         s=[50.0], alpha=.8)

        if handles is None:
            handles, labels = plot.get_legend_handles_labels()

    plt.figlegend(labels=labels,
                  loc="lower center", handles=handles,
                  scatterpoints=1,
                  ncol=ceil(len(options["methods"])*2/4), fontsize=10,
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
        plt.savefig("{}.{}".format(options["out"], options["format"]),
                    format=options["format"],
                    dpi=options["dpi"],
                    transparent=options["opaque"])

if __name__ == "__main__":
    main()
