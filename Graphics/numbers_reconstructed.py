import matplotlib.pyplot as plt
from math import ceil
from collections import OrderedDict
import matplotlib.lines as mlines
import matplotlib.colors
import matplotlib.cm as cm
import csv
import argparse
from itertools import zip_longest
import matplotlib.patches as mpatches
from utils import parse_configuration
import os
import re


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def clamp(x):
    return max(0, min(x, 255))


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--title", required=True)
    parser.add_argument("--log", action="store_true", default=False)
    parser.add_argument("--out", required=True)
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    data = OrderedDict()

    options = parse_configuration(args.configuration)

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(options["colourmap"]["name"])

    for label in options["methods"]:
        for aligner in options["divisions"]:
            if options["methods"][label][aligner] is not None:
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
            else:
                data[(label, aligner)] = [-1000] * 3

    # print(*data.items(), sep="\n")

    divisions = sorted(options["divisions"].keys())

    figure, axes = plt.subplots(nrows=3,
                                ncols=1,
                                dpi=300, figsize=(8, 6))
    figure.suptitle(args.title)

    newticks = ["Recovered genes", "Missed genes", "Fused genes"]
    for pos, ax in enumerate(axes):
        ax.set_ylim(0, 2)
        max_x = max(data[_][pos] for _ in data) + 500
        min_x = max(0, min(data[_][pos] for _ in data if data[_][pos] > 0) - 500)
        ax.set_xlim(min_x, max_x)
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

    handles = []
    labels = []
    for index, tup in enumerate(data.keys()):
        method, division = tup
        if options["colourmap"]["use"] is False:
            colour = options["methods"][method]["colour"]
            matched = re.match("\(([0-9]*), ([0-9]*), ([0-9]*)\)$", colour)
            if matched:
                colour = "#{0:02x}{1:02x}{2:02x}".format(clamp(int(matched.groups()[0])),
                                                         clamp(int(matched.groups()[1])),
                                                         clamp(int(matched.groups()[2])))
        elif options["methods"][method]["colour"] in ("black", "k"):
            colour = "black"
        else:
            colour = color_map(color_normalizer(options["methods"][division]["index"]))

        # color = colors[int(index / 2)]
        marker = options["divisions"][division]["marker"]
        cat = "{} ({})".format(method, division)
        labels.append(cat)
        # handles.append(handle)
        for pos, point in enumerate(data[tup]):
            handle = axes[pos].scatter(point,
                                        1,
                                        alpha=1,
                                        label=cat,
                                        color=colour,
                                        marker=marker,
                                        edgecolor="k",
                                        s=100)
            handle.get_sketch_params()
            if pos == 0:
                handles.append(handle)
        handle = mlines.Line2D([], [], markersize=5, color=colour, marker=marker, label=cat, alpha=0.6)

    div_labels = []
    for division in options["divisions"]:
        faux_line = mlines.Line2D([], [], color="white",
                                  marker=options["divisions"][division]["marker"],
                                  markersize=14,
                                  markerfacecolor="black")
        div_labels.append((faux_line, division))

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

    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="lower center",
                  scatterpoints=1,
                  ncol=ceil(len(options["methods"])*2/4), fontsize=10,
                  framealpha=0.5)
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
