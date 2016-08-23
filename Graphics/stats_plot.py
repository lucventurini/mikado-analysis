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
from math import ceil, floor
from itertools import zip_longest

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
    parser.add_argument("--tophat", nargs="+", required=True,
                        help="For each invocation, specify the original and the filtered stats file from Compare.")
    parser.add_argument("--star", nargs="+", required=True,
                        help="For each invocation, specify the original and the filtered stats file from Compare.")
    parser.add_argument("--labels", nargs="+", required=True,
                        help="Labels to use. They must be in the same number and order of the star/tophat files.")
    parser.add_argument("--out", default=None, required=True)
    parser.add_argument("--format", default="svg", choices=["svg",
                                                            "pdf",
                                                            "png"])
    parser.add_argument("-c", "--colours", "--colors", dest="colours",
                        default=None, nargs="+",
                        help="Colours to be used. Defaults to use a colourmap")
    parser.add_argument("-cm", "--colourmap", default="gist_rainbow",
                        help="Colourmap to be used.")
    parser.add_argument("--dpi", default=600, type=int)
    parser.add_argument("--title", default="Mikado stats")
    args = parser.parse_args()

    if (len(args.tophat) != len(args.star) or len(args.labels) != int(len(args.star) / 2) or
            len(args.star) % 2 != 0):
        print(
            "Error, labels and stats files are not the same number \
            (TH: {}, STAR: {}, Labels: {})\n\
             TH: {}\n\
             STAR: {}\n\
            Labels: {}".format(
                len(args.tophat), len(args.star), len(args.labels),
            args.tophat, args.star, args.labels))
        parser.print_help()
        sys.exit(1)

    if args.colours is not None:
        if len(args.colours) != len(args.labels) or len(set(args.colours)) != len(args.colours):
            print("Error, invalid number of unique colors specified")
            parser.print_help()
            sys.exit(1)

    stats = OrderedDict()
    figure, axes = plt.subplots(nrows=2,
                                ncols=3,
                                dpi=args.dpi, figsize=(8, 6))
    figure.suptitle("${}$".format(args.title),
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

    name_ar = np.array([["Base", "Exon", "Intron"],
                        ["Intron chain", "Transcript", "Gene"]])

    color_normalizer = matplotlib.colors.Normalize(0, len(args.labels))
    color_map = cm.get_cmap(args.colourmap)
    # mapper = cm.ScalarMappable(colors, "PuOr")

    for xrow in (0, 1):
        for yrow in (0, 1, 2):
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

    print(*args.star, sep="\n")
    
    # We presume we have FIRST all original, THEN all filtered.
    args.star = list(zip(args.star[:len(args.labels)], args.star[len(args.labels):]))
    args.tophat = list(zip(args.tophat[:len(args.labels)], args.tophat[len(args.labels):]))

    for label, name in zip(args.labels, args.star):
        print(label, name)
        orig, filtered = name
        # orig = "{}-compare.stats".format(name)
        # filtered = "{}-filtered_compare.stats".format(name)

        orig_lines = [line.rstrip() for line in open(orig)]
        filtered_lines = [line.rstrip() for line in open(filtered)]
        # In the stats we have precision as second and sensitivity as first,
        # we have to invert
        for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            precision = float(orig_lines[line_index].split(":")[1].split()[1])
            recall = float(filtered_lines[line_index].split(":")[1].split()[0])
            stats[
                # Name of the statistic:Base, Exon, etc
                list(stats.keys())[index]][
                b"STAR"].append((precision, recall))

    for label, name in zip(args.labels, args.tophat):
        print(label, name)
        orig, filtered = name
        # orig = "{}-compare.stats".format(name)
        # filtered = "{}-filtered_compare.stats".format(name)

        orig_lines = [line.rstrip() for line in open(orig)]
        filtered_lines = [line.rstrip() for line in open(filtered)]
        # In the stats we have precision as second and sensitivity as first,
        # we have to invert
        for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            precision = float(orig_lines[line_index].split(":")[1].split()[1])
            recall = float(filtered_lines[line_index].split(":")[1].split()[0])
            stats[
                # Name of the statistic:Base, Exon, etc
                list(stats.keys())[index]][
                b"TopHat"].append((precision, recall))

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

        for index, vals in enumerate(zip(Xtop, Ytop, args.labels)):
            x, y, label = vals
            if args.colours is not None:
                colour = args.colours[index]
            else:
                colour = color_map(color_normalizer(index))
            plot.scatter(x, y, label="{0} (TopHat)".format(label), c=colour,
                         marker="^", edgecolor="k",
                         s=[50.0], alpha=.8)

        for index, vals in enumerate(zip(Xstar, Ystar, args.labels)):
            x, y, label = vals
            if args.colours is not None:
                colour = args.colours[index]
            else:
                colour = color_map(color_normalizer(index))
            plot.scatter(x, y, label="{0} (STAR)".format(label), c=colour,
                         marker="o",
                         s=[50.0], alpha=.8)

        if handles is None:
            handles, labels = plot.get_legend_handles_labels()

    plt.figlegend(labels=labels,
                  loc="lower center", handles=handles,
                  scatterpoints=1,
                  ncol=3, fontsize=10)
    # Necessary to pad the superior title
    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.1,  # Left
                           0.2,  # Bottom
                           0.85,  # Right
                           0.9])  # Top
    if args.out is None:
        plt.show()
    else:
        plt.savefig(args.out, format=args.format, dpi=args.dpi)

if __name__ == "__main__":
    main()
