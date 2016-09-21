import matplotlib.pyplot as plt
from collections import OrderedDict
import matplotlib.lines as mlines
import csv
import argparse
from itertools import zip_longest
import re
from utils import parse_configuration
import os


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--title", required=True)
    parser.add_argument("--log", action="store_true", default=False)
    parser.add_argument("--out", required=True)
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    data = OrderedDict()

    options = parse_configuration(args)

    for label in options["methods"]:
        for aligner in ("STAR", "TopHat"):
            data[(label, aligner)] = [set(), set(), set()]
            orig_refmap = "{}.refmap".format(
                re.sub(".stats$", "", options["methods"][label][aligner][0]))
            with open(orig_refmap) as refmap:
                for row in csv.DictReader(refmap, delimiter="\t"):
                    if row["best_ccode"] in ("=", "_"):
                        data[(label, aligner)][0].add(row["ref_gene"])
                    elif row["best_ccode"][0] == "f":
                        data[(label, aligner)][2].add(row["ref_gene"])
                    # elif row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                    #     data[(label, aligner)][1].add(row["ref_gene"])
            filtered_refmap = "{}.refmap".format(
                re.sub(".stats$", "", options["methods"][label][aligner][1]))
            with open(filtered_refmap) as refmap:
                for row in csv.DictReader(refmap, delimiter="\t"):
                    if row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                        data[(label, aligner)][1].add(row["ref_gene"])
            for num in range(3):
                data[(label, aligner)][num] = len(data[(label, aligner)][num])

    # print(*data.items(), sep="\n")

    figure, axes = plt.subplots(nrows=3,
                                ncols=1,
                                dpi=300, figsize=(8, 6))
    figure.suptitle(args.title)

    newticks = ["Recovered genes", "Missed genes", "Fused genes"]
    for pos, ax in enumerate(axes):
        ax.set_ylim(0, 2)
        max_x = max(data[_][pos] for _ in data) * 1.1
        ax.set_xlim(0, max_x)
        ax.plot((1, max_x), (1, 1), 'k-')
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.set_yticks([1])
        ax.set_yticklabels([newticks[pos]], fontsize=15)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", left="off", right="off")
        ax.tick_params(axis="x", top="off")
        if pos == 2:
            ax.set_xlabel("Number of genes", fontsize=12)

    # Set the axis to log if necessary
    if args.log is True:
        plt.xscale("log")

    #
    # plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (2, 2), 'k-')
    # plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (3, 3), 'k-')

    shape = ["o", "^"]

    handles = []
    labels = []
    for index, tup in enumerate(data.keys()):

        method, aligner = tup
        color = options["methods"][method]["colour"]
        # color = colors[int(index / 2)]
        marker = shape[index % 2]
        cat = "{} ({})".format(method, aligner)
        labels.append(cat)
        # handles.append(handle)
        for pos, point in enumerate(data[tup]):
            handle = axes[pos].scatter(point,
                                        1,
                                        alpha=1,
                                        label=cat,
                                        color=color,
                                        marker=marker,
                                        edgecolor="k",
                                        s=100)
            handle.get_sketch_params()
            if pos == 0:
                handles.append(handle)
        handle = mlines.Line2D([], [], markersize=5, color=color, marker=marker, label=cat, alpha=0.6)

    # handles, labels = plot.get_legend_handles_labels()
    plt.figlegend(labels=labels,
                  framealpha=0.3,
                  loc=(0.1, 0.09), handles=handles,
                  scatterpoints=1,
                  ncol=3, fontsize=10, markerscale=0.6)
    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.1,  # Left
                           0.25,  # Bottom
                           1,  # Right
                           0.9])  # Top
    out = "{}.{}".format(os.path.splitext(args.out)[0], options["format"])
    plt.savefig(out,
                format=options["format"],
                transparent=True)

main()
