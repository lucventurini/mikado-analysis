#!/usr/bin/env python3

import argparse
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
    parser.add_argument("-cm", "--colour-map", dest="colour_map",
                        default=None,
                        help="Colour map to use. Default: use user-specified colours.")
    parser.add_argument("-cs", "--colour-map-size", dest="colour_map_size",
                        default=0,
                        type=int,
                        help="Colour map to use. Default: use user-specified colours.")
    parser.add_argument("--title", default="Mikado stats for transcript level")
    args = parser.parse_args()

    species = yaml.load(args.species)

    ncols = ceil(len(species) / 2)

    figure, axes = plt.subplots(
        nrows=2,
        ncols=ncols,
        dpi=args.dpi,
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

    name_ar = np.array(list(grouper(species.keys(), ceil(len(species)/2), None)))

    # Dictionary to indicate which line should be taken given the level
    line_correspondence = {"base": 5,
                           "exon": 7,
                           "intron": 8,
                           "intron_chain": 9,
                           "transcript": 12,
                           "gene": 15}

    markers = dict()
    methods = dict()

    best_marker = "x"

    for xrow in (0, 1):
        for yrow in range(ncols):
            try:
                key = name_ar[xrow, yrow]
            except IndexError:
                raise IndexError(name_ar)
            if key is None:
                continue
            plot = axes[xrow, yrow]
            # plot.grid(True, linestyle='dotted')
            # plot.set(adjustable="box-forced", aspect="equal")
            plot.set_title("${}$".format(species[key]["name"]), fontsize=10)

            # plot.set_xlabel("Precision", fontsize=10)
            # plot.set_ylabel("Recall", fontsize=10)

            stats = OrderedDict()

            with open(species[key]["configuration"]) as configuration:
                options = parse_configuration(configuration, prefix=species[key]["folder"])
                if args.colour_map is not None:
                    options["colourmap"]["use"] = True
                    options["colourmap"]["name"] = args.colour_map
                else:
                    options["colourmap"]["use"] = False

                for division in options["divisions"]:
                    markers[division] = options["divisions"][division]["marker"]

                if options["colourmap"]["use"] is True:
                    cm_length = max(args.colour_map_size, len(options["methods"]))
                    print(cm_length)
                    color_normalizer = matplotlib.colors.Normalize(vmin=0, vmax=cm_length)
                    color_map = cm.get_cmap(options["colourmap"]["name"])

                for division in options["divisions"]:
                    stats[division.encode()] = []

                for method in options["methods"]:
                    for aligner in options["divisions"]:
                        orig, filtered = options["methods"][method][aligner]
                        orig_lines = [line.rstrip() for line in open(orig)]
                        filtered_lines = [line.rstrip() for line in open(filtered)]
                        # for index, line_index in enumerate([5, 7, 8, 9, 11, 12, 14, 15]):
                        for index, line_index in enumerate([line_correspondence[args.level]]):
                            precision = float(orig_lines[line_index].split(":")[1].split()[1])
                            recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                            try:
                                f1 = hmean(np.array([precision, recall]))
                            except TypeError as exc:
                                raise TypeError("\n".join([str(_) for _ in [(precision, type(precision)),
                                                                            (recall, type(recall)),
                                                                            exc]]))
                            stats[aligner.encode()].append((precision, recall, f1))
                divisions = sorted(options["divisions"].keys())
                handles = None
                ys = [np.array([_[0] for _ in stats[division.encode()]]) for division in divisions]
                xs = [np.array([_[1] for _ in stats[division.encode()]]) for division in divisions]
                x_minimum = max(0,
                                floor(min(_.min() for _ in xs)) - 5)
                y_minimum = max(0,
                                floor(min(_.min() for _ in ys)) - 5)

                x_maximum = min(100,
                                ceil(max(_.max() for _ in xs)) + 5)
                y_maximum = min(100,
                                ceil(max(_.max() for _ in ys)) + 5)
                plotf1curves(plot, fstepsize=ceil(min(x_maximum - x_minimum, y_maximum - y_minimum) / 10))
                best_f1 = (-1, [])

                for enumerated, division in enumerate(divisions):
                    for index, vals in enumerate(zip(xs[enumerated], ys[enumerated], options["methods"].keys())):
                        x, y, label = vals

                        f1 = calc_f1(x, y)
                        if best_f1[0] < f1:
                            best_f1 = (f1, [(x, y)])
                        elif best_f1[0] == f1:
                            best_f1[1].append((x, y))

                        if label in ("Mikado permissive", "Mikado stringent"):
                            method_name = "Mikado"
                        else:
                            method_name = label
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
                        methods[method_name] = colour
                        plot.scatter(x, y,
                                     label=label,
                                     # label="{0} ({1})".format(label, division),
                                     c=colour, marker=options["divisions"][division]["marker"],
                                     edgecolor="k", s=[50.0], alpha=1)
                __axes = plot.axes
                __axes.set_xlim(x_minimum, x_maximum)
                __axes.set_ylim(y_minimum, y_maximum)
                # __axes.set_aspect("equal")
                plot.tick_params(axis='both', which='major', labelsize=8)
                # Annotate the best F1
                for best in best_f1[1]:
                    plot.scatter(best[0], best[1],
                                 label="Best F1",
                                 marker=best_marker, s=[20], c="k")

    div_labels = []

    f1_line = mlines.Line2D([], [], color="gray", linestyle="--")
    div_labels.append((f1_line, "F1 contours"))

    for division in markers:
        faux_line = mlines.Line2D([], [], color="white",
                                  marker=markers[division],
                                  markersize=14,
                                  markerfacecolor="black")
        div_labels.append((faux_line, division))

    best_marker_line = mlines.Line2D([], [], color="white",
                                marker=best_marker, markersize=6,
                                markerfacecolor="black", markeredgecolor="black")
    div_labels.append((best_marker_line, "Best F1"))

    for method in methods:
        colour = methods[method]
        patch = mpatches.Patch(facecolor=colour, linewidth=1, edgecolor="k")
        div_labels.append((patch, method))

    plt.figlegend(handles=[_[0] for _ in div_labels],
                  labels=[_[1] for _ in div_labels],
                  loc="lower center",
                  scatterpoints=1,
                  ncol=ceil((len(methods) + len(markers)) / 4),
                  fontsize=10,
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