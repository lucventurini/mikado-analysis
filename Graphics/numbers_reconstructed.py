import matplotlib.pyplot as plt
from math import ceil
from operator import itemgetter
from collections import OrderedDict, defaultdict
import matplotlib.lines as mlines
import matplotlib.colors
import matplotlib.cm as cm
import csv
import argparse
from itertools import zip_longest
import matplotlib.patches as mpatches
from utils import parse_configuration, parse_refmaps
import os
import re
import sys
import numpy as np
import pandas as pd
import multiprocessing
# from brokenaxes import brokenaxes


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
    parser.add_argument("--title", required=False, default="")
    parser.add_argument("--log", action="store_true", default=False)
    parser.add_argument("--out", required=False, default=None, help="Output file. If unspecified, the script will exit after printing the numbers.")
    parser.add_argument("--genes", default=False, action="store_true",
                        help="Flag. If switched on, the gene level will be used instead of the transcript level.")
    parser.add_argument("-p", "--procs", default=multiprocessing.cpu_count(), type=int)
    parser.add_argument("--transcripts", action="store_true", default=False)
    parser.add_argument("--division", action="store_true", default=False)
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    data = OrderedDict()

    options = parse_configuration(args.configuration)

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(options["colourmap"]["name"])

    header = ["Method", "Division", "Fully", "Missed", "Fused"]
    print(*header, sep="\t")

    pool = multiprocessing.Pool(processes=args.procs)

    for label in options["methods"]:
        for aligner in options["divisions"]:
            if options["methods"][label][aligner] is not None:
                data[(label, aligner)] = [set(), set(), set()]
                orig_refmap = "{}.refmap".format(
                    re.sub(".stats$", "", options["methods"][label][aligner][0]))
                with open(orig_refmap) as refmap:
                    for row in csv.DictReader(refmap, delimiter="\t"):
                        if args.genes is True:
                            if row["best_ccode"] in ("=", "_"):
                                data[(label, aligner)][0].add(row["ref_gene"])
                            elif row["best_ccode"][0] == "f":
                                data[(label, aligner)][2].add(row["ref_gene"])
                            # elif row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                            #     data[(label, aligner)][1].add(row["ref_gene"])
                        else:
                            if row["ccode"] in ("=", "_"):
                                data[(label, aligner)][0].add(row["ref_id"])
                            elif row["ccode"][0] == "f":
                                data[(label, aligner)][2].add(row["ref_id"])
                filtered_refmap = "{}.refmap".format(
                    re.sub(".stats$", "", options["methods"][label][aligner][1]))
                with open(filtered_refmap) as refmap:
                    for row in csv.DictReader(refmap, delimiter="\t"):
                        if args.genes is True:
                            if row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                                data[(label, aligner)][1].add(row["ref_gene"])
                        else:
                            if row["ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                                data[(label, aligner)][1].add(row["ref_id"])
                for num in range(3):
                    data[(label, aligner)][num] = len(data[(label, aligner)][num])
                orig_stats, filtered_stats = options["methods"][label][aligner][:2]
            else:
                orig_stats, filtered_stats = None, None
            data[(label, aligner)] = pool.apply_async(parse_refmaps, (orig_stats, filtered_stats, args.transcripts))

    for label in options["methods"]:
        for aligner in options["divisions"]:
            data[(label, aligner)] = data[(label, aligner)].get()
            print(label, aligner, *data[(label, aligner)], sep="\t")

    # print(*data.items(), sep="\n")

    # Now print out the table

    if args.out is None:
        sys.exit(0)

    # divisions = sorted(options["divisions"].keys())

    if args.division is True:
        figsize = (14, 12)
    else:
        figsize = (6, 8)

    figure, axes = plt.subplots(nrows=1,
                                ncols=3,
                                dpi=300, figsize=figsize)
    figure.suptitle(args.title)

    handles = []
    labels = []

    factor = 10

    for pos, ax in enumerate(axes):
        ax.set_ylim(0, 2)
        max_x = max(data[_][pos] for _ in data) + 500
        min_x = max(0, min(data[_][pos] for _ in data if data[_][pos] > 0) - 500)
        ax.set_ylim(min_x, max_x)
        if args.division is False:
            ax.set_xlim(0, 2.5)
        else:
            ax.set_xlim(-2, len(options["divisions"]) * factor)
        # ax.plot((1, max_x), (1, 1), 'k-')
        ax.tick_params(axis='both', which='major', labelsize=10)
        # ax.set_xticklabels([newticks[pos]], fontsize=15)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["bottom"].set_visible(True)
        ax.tick_params(axis="y", left="on", right="off")
        if args.division is False:
            ax.tick_params(axis="x", top="off", bottom="off")
        else:
            ax.tick_params(axis="x", top="off", bottom="on")
        if pos == 0:
            if args.transcripts is True:
                ax.set_ylabel("Number of transcripts", fontsize=16)
            else:
                ax.set_ylabel("Number of genes", fontsize=16)
            ax.set_xlabel("Reconstructed\ngenes", fontsize=16)
        elif pos == 1:
            ax.set_xlabel("Missed\ngenes", fontsize=16)
        else:
            ax.set_xlabel("Fused\ngenes", fontsize=16)

        if args.division is True:
            ax.set_xticks([factor * (_ + 1 / 4) for _ in range(len(options["divisions"]))])
            ax.set_xticklabels(options["divisions"], rotation=45, fontsize=14)

        if args.division is False:
            points = []
        else:
            points = defaultdict(list)

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

            marker = options["divisions"][division]["marker"]
            cat = "{} ({})".format(method, division)
            labels.append(cat)

            point = data[tup][pos]
            if args.division is False:
                points.append([point, cat, colour, marker])
            else:
                points[division].append([point, cat, colour, marker])

        if args.division is False:
            points = sorted(points, key=itemgetter(0))
            for index, point in enumerate(points):
                point, cat, colour, marker = point
                x_coord = 0.5 + (index % 4 / 2)
                handle = axes[pos].scatter(x_coord, point,
                                           alpha=1,
                                           label=cat,
                                           color=colour,
                                           marker=marker,
                                           edgecolor="k",
                                           s=150)
                handle.get_sketch_params()
                if pos == 0:
                    handles.append(handle)
                handle = mlines.Line2D([], [], markersize=5, color=colour, marker=marker, label=cat, alpha=0.6)
        else:

            max_last_xcoord = None

            for dindex, division in enumerate(options["divisions"].keys()):
                max_xcoord = -100
                min_xcoord = 100000
                dpoints = sorted(points[division], key=itemgetter(0))
                for index, point in enumerate(dpoints):
                    point, cat, colour, marker = point
                    if index % 3 == 0:
                        x_coord = factor * dindex
                    elif index % 3 == 2:
                        x_coord = factor * (dindex + 1 / 2)
                    else:
                        x_coord = factor * (dindex + 1 / 4)

                    max_xcoord = max(x_coord, max_xcoord)
                    min_xcoord = min(x_coord, min_xcoord)
                    # print(x_coord, point, cat)
                    handle = axes[pos].scatter(x_coord, point,
                                               alpha=1,
                                               label=cat,
                                               color=colour,
                                               marker=marker,
                                               edgecolor="k",
                                               s=150)
                    handle.get_sketch_params()
                    if pos == 0:
                        handles.append(handle)
                    handle = mlines.Line2D([], [], markersize=5, color=colour, marker=marker, label=cat, alpha=0.6)
                if max_last_xcoord is not None and min_xcoord < max_last_xcoord:
                    raise ValueError("Overlapping X values")
                max_last_xcoord = max_xcoord

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
                                        s=150)
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
        elif options["methods"][method]["colour"] in ("black", "k"):
            colour = "black"
        else:
            colour = color_map(color_normalizer(options["methods"][method]["index"]))

        patch = mpatches.Patch(facecolor=colour, linewidth=1, edgecolor="k")
        div_labels.append((patch, method))

    if args.division is True:
        fs = 16
    else:
        fs = 10

    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="lower center",
                  scatterpoints=1,
                  ncol=min(ceil(len(options["methods"])*2/4), 3), fontsize=fs,
                  framealpha=0.5)
    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.05,  # Left
                           0.15,  # Bottom
                           0.95,  # Right
                           0.95])  # Top
    out = "{}.{}".format(os.path.splitext(args.out)[0], options["format"])
    plt.savefig(out,
                format=options["format"],
                transparent=True)

main()
