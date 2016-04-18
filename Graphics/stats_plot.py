#!/usr/bin/env python3

import argparse
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
from collections import OrderedDict
import matplotlib.ticker as ticker


__doc__ = """Script to automate the plots for the Mikado compare statistics"""



def split_comma(string):
    return string.split(",")



def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--labels", type=split_comma, required=True)
    parser.add_argument("--stats", type=split_comma, required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--format", default="svg", choices=["svg",
                                                            "pdf",
                                                            "png"])
    parser.add_argument("--dpi", default=600, type=int)
    parser.add_argument("--title", default="Mikado stats")
    args = parser.parse_args()

    if len(args.labels) != len(args.stats):
        print("Error, labels and stats files are not the same number")
        parser.print_help()
        sys.exit(1)

    stats = OrderedDict()
    figure, axes = plt.subplots(nrows=2,
                                ncols=3,
                                dpi=args.dpi, figsize=(10, 6))
    figure.suptitle(args.title, fontsize=15)
    # plt.setp(axes, xticks=[0.1, 0.5, 0.9], xticklabels=['a', 'b', 'c'],
    #          yticks=[1, 2, 3])

    name_ar = np.array([["Base", "Exon", "Intron"],
                        ["Intron chain", "Transcript", "Gene"]])

    markers = ["o", "v", "8", "s", "p", "H", "+", "x"]

    for position in name_ar:
        print(position)

    color_normalizer = matplotlib.colors.Normalize(0, len(args.labels))
    color_map = cm.get_cmap("Spectral")
    # mapper = cm.ScalarMappable(colors, "PuOr")

    for xrow in (0, 1):
        for yrow in (0, 1, 2):
            key = name_ar[xrow, yrow]
            plot = axes[xrow, yrow]
            plot.grid(True)
            plot.set_title("{} level".format(key), fontsize=12)
            plot.axis("scaled")
            plot.set_xlim(0, 100)
            plot.set_ylim(0, 100)

            plot.plot(plot.get_xlim(), plot.get_ylim(), ls="--", c=".3")

            __axes = plot.axes
            __axes.xaxis.set_major_locator(ticker.MultipleLocator(30))
            __axes.xaxis.set_minor_locator(ticker.MultipleLocator(5))
            __axes.yaxis.set_major_locator(ticker.MultipleLocator(30))
            __axes.yaxis.set_minor_locator(ticker.MultipleLocator(5))
            plot.set_xlabel("Precision", fontsize=6)
            plot.set_ylabel("Recall", fontsize=6)
            stats[key] = dict()
            stats[key][b"plot"] = plot
            stats[key][b"stats"] = []

    for name in args.stats:
        lines = [line.rstrip() for line in open(name)]
        # In the stats we have precision as second and sensitivity as first,
        # we have to invert
        for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            stats[list(stats.keys())[index]][b"stats"].append(
                (list(reversed([float(_) for _ in lines[line_index].split(":")[1].split()[:2]])))
            )

    handles, labels = None, None
    for stat in stats.keys():
        values = stats[stat][b"stats"]
        plot = stats[stat][b"plot"]

        X = np.array([_[0] for _ in values])
        Y = np.array([_[1] for _ in values])
        for index, vals in enumerate(zip(X, Y, args.labels)):
            x, y, label = vals
            color = color_map(color_normalizer(index))
            plot.scatter(x, y, label=label, c=color,
                         marker=markers[index],
                         s=[100.0], alpha=.6)

        # plot.scatter(X, Y)
        # plot.legend(loc="lower center", fontsize=5, ncol=2)
        if handles is None:
            handles, labels = plot.get_legend_handles_labels()
            # print(handles)
            # print(labels)

    # print(handles)
    plt.figlegend(labels=labels,
                  loc="lower center", handles=handles,
                  scatterpoints=1,
                  ncol=3, fontsize=10)
    # Necessary to pad the superior title
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    if args.out is None:
        plt.show()
    else:
        plt.savefig(args.out, format=args.format, dpi=600)

main()