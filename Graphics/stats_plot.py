#!/usr/bin/env python3

import argparse
import sys
import matplotlib.lines
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
from collections import OrderedDict
import matplotlib.ticker as ticker
from math import ceil


__doc__ = """Script to automate the plots for the Mikado compare statistics"""


def split_comma(string):
    return string.split(",")


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--tophat", type=split_comma, required=True)
    parser.add_argument("--star", type=split_comma, required=True)
    parser.add_argument("--labels", type=split_comma, required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--format", default="svg", choices=["svg",
                                                            "pdf",
                                                            "png"])
    parser.add_argument("--dpi", default=600, type=int)
    parser.add_argument("--title", default="Mikado stats")
    args = parser.parse_args()

    if len(args.labels) != len(args.tophat) != len(args.star):
        print("Error, labels and stats files are not the same number")
        parser.print_help()
        sys.exit(1)

    stats = OrderedDict()
    figure, axes = plt.subplots(nrows=2,
                                ncols=3,
                                dpi=args.dpi, figsize=(10, 6))
    figure.suptitle(args.title, fontsize=15)

    Xaxis = matplotlib.patches.FancyArrow(0.06, 0.16, 0.9, 0,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")
    Yaxis = matplotlib.patches.FancyArrow(0.06, 0.16, 0, 0.75,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")
    figure.lines.extend([Xaxis, Yaxis])

    figure.text(0.9, 0.13, "$Precision$", ha="center", fontsize=12)
    figure.text(0.04, 0.75, "$Recall$", va="center", fontsize=12, rotation="vertical")


    # plt.setp(axes, xticks=[0.1, 0.5, 0.9], xticklabels=['a', 'b', 'c'],
    #          yticks=[1, 2, 3])

    name_ar = np.array([["Base", "Exon", "Intron"],
                        ["Intron chain", "Transcript", "Gene"]])

    color_normalizer = matplotlib.colors.Normalize(0, len(args.labels))
    color_map = cm.get_cmap("gist_rainbow")
    # mapper = cm.ScalarMappable(colors, "PuOr")

    for xrow in (0, 1):
        for yrow in (0, 1, 2):
            key = name_ar[xrow, yrow]
            plot = axes[xrow, yrow]
            plot.grid(True)
            plot.set_title("{} level".format(key), fontsize=12)

            # plot.set_xlabel("Precision", fontsize=10)
            # plot.set_ylabel("Recall", fontsize=10)
            stats[key] = dict()
            stats[key][b"plot"] = plot
            stats[key][b"STAR"] = []
            stats[key][b"TopHat"] = []

    for name in args.star:
        lines = [line.rstrip() for line in open(name)]
        # In the stats we have precision as second and sensitivity as first,
        # we have to invert
        for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            stats[
                # Name of the statistic:Base, Exon, etc
                list(stats.keys())[index]][
                b"STAR"].append(
                (list(reversed([float(_) for _ in lines[line_index].split(":")[1].split()[:2]])))
            )

    for name in args.tophat:
        lines = [line.rstrip() for line in open(name)]
        # In the stats we have precision as second and sensitivity as first,
        # we have to invert
        for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            stats[
                # Name of the statistic:Base, Exon, etc
                list(stats.keys())[index]][
                b"TopHat"].append(
                # Append a tuple of values, Precision and Recall
                (list(reversed([float(_) for _ in lines[line_index].split(":")[1].split()[:2]])))
            )

    handles, labels = None, None

    for stat in stats.keys():
        star_values = stats[stat][b"STAR"]
        th_values = stats[stat][b"TopHat"]
        plot = stats[stat][b"plot"]

        Xtop = np.array([_[0] for _ in th_values])
        Ytop = np.array([_[1] for _ in th_values])

        Xstar = np.array([_[0] for _ in star_values])
        Ystar = np.array([_[1] for _ in star_values])

        plot.axis("scaled")
        # Select a suitable maximum

        maximum = 10 * ceil(
            min(100,
                max(Xtop.max(), Xstar.max(), Ystar.max(), Ytop.max())
                ) / 10.0 )

        plot.set_xlim(0, maximum)
        plot.set_ylim(0, maximum)
        __axes = plot.axes
        __axes.xaxis.set_major_locator(ticker.MultipleLocator(ceil(maximum / 20) * 5))
        __axes.xaxis.set_minor_locator(ticker.MultipleLocator(ceil(maximum / 100) * 5))
        __axes.yaxis.set_major_locator(ticker.MultipleLocator(ceil(maximum / 20) * 5))
        __axes.yaxis.set_minor_locator(ticker.MultipleLocator(ceil(maximum / 100) * 5))

        plot.plot(plot.get_xlim(), plot.get_ylim(), ls="--", c=".3")


        for index, vals in enumerate(zip(Xtop, Ytop, args.labels)):
            x, y, label = vals
            colour = color_map(color_normalizer(index))
            plot.scatter(x, y, label="{0} (TopHat)".format(label), c=colour,
                         marker="^", edgecolor="k",
                         s=[100.0], alpha=.8)

        for index, vals in enumerate(zip(Xstar, Ystar, args.labels)):
            x, y, label = vals
            colour = color_map(color_normalizer(index))
            plot.scatter(x, y, label="{0} (STAR)".format(label), c=colour,
                         marker="o",
                         s=[100.0], alpha=.8)

        if handles is None:
            handles, labels = plot.get_legend_handles_labels()

    plt.figlegend(labels=labels,
                  loc="lower center", handles=handles,
                  scatterpoints=1,
                  ncol=4, fontsize=10)
    # Necessary to pad the superior title
    plt.tight_layout(rect=[0, 0.15, 1, 0.95])
    if args.out is None:
        plt.show()
    else:
        plt.savefig(args.out, format=args.format, dpi=600)

main()