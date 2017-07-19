#!/usr/bin/env python3

import sys
import argparse
import csv
import functools
import re
import pandas
from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches
from math import floor
import seaborn
from utils import parse_configuration
pandas.options.mode.chained_assignment = None

def sort_values(key, dictionary):
    return dictionary[key]["TPM"]


def generate_plot(dataframe, args, options, nrows=2, ncols=5):

    plt.style.context("ggplot")
    dataframe = dataframe.sort_values("TPM", ascending=True)

    order = {0: [dataframe[(dataframe.TPM <= 0.01)],
                      "$\leq$ $0.01$"],
             1: [dataframe[(dataframe.TPM > 0.01) & (dataframe.TPM <= 2)],
                      "$1$ $-$ $2$"],
             2: [dataframe[(dataframe.TPM > 2) & (dataframe.TPM <= 10)],
                      "$2$ $-$ $10$"],
             3: [dataframe[(dataframe.TPM <= 100) & (dataframe.TPM > 10)],
                      "$10$ $-$ $100$"],
             4: [dataframe[(dataframe.TPM > 100)],
                      "$\geq$ 100"]}

    color_map = cm.get_cmap(options["colourmap"]["name"])
    color_normalizer = matplotlib.colors.Normalize(2, 15)

    colors = {1: "white", 2: "white"}

    for num in range(3, 8 * 2 + 1):
        colors[num] = color_map(color_normalizer(num))

    # plot = plt.plot()

    figure, axes=plt.subplots(nrows=nrows,
                              ncols=ncols,
                              figsize=(16, 10))

    methods = dataframe.columns[2:]
    current = 0

    figure.suptitle(" ".join(["${}$".format(re.sub("%", "\%", _)) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")

    Xaxis = matplotlib.patches.FancyArrow(0.05, 0.1, 0.89, 0,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")

    Yaxis = matplotlib.patches.FancyArrow(0.05, 0.1, 0, 0.8,
                                          width=0.0005,
                                          length_includes_head=True,
                                          transform=figure.transFigure, figure=figure,
                                          color="k")
    figure.lines.extend([Xaxis, Yaxis])

    figure.text(0.92, 0.12, "$Expression$ $(TPM)$", ha="center", fontsize=15)
    figure.text(0.02,   # X location
                0.6,  # Y location
                "$Transcripts$ $per$ $category$ $(\%)$",  # Title
                va="center",
                fontsize=15,
                rotation="vertical")

    legend_handles = []
    current = 0

    # Structure of the dataset:
    # <TID> <Expr> <Method1> <Method2> ... <MethodN>

    for row, aligner in enumerate(sorted(options["divisions"])):
        for col in range(5):
            # if current == 6:
            #     axes[row, col].xaxis.set_visible(False)
            #     axes[row, col].yaxis.set_visible(False)
            #     continue
            fraction = order[col][0]
            if nrows > 1:
                plot = axes[row, col]
            else:
                plot = axes[col]
            plot.set_title("{} $TPM$ (${}$)".format(order[col][1], aligner))
            plot.set_ylim(0, 100)
            _axes = plot.axes
            _axes.set_xticks(numpy.arange(0, 5, 1))
            curr_labels = []

            for tick in plot.get_xticklabels():
                tick.set_rotation(60)
            values_array = []

            __vals = fraction[[_ for n, _ in enumerate(fraction.columns) if n > 1]]

            __vals.dropna(inplace=True)
            __vals = __vals[__vals.apply(lambda x: min(x) == 6, 1)]
            found_in_all = len(__vals)

            for method in methods:
                if aligner not in method:
                    continue
                # Get for each method the number of transcripts
                # Corresponding to each ccode
                curr_labels.append("\n".join([
                    "${}$".format(_) for _ in
                    re.sub(" \({}\)".format(aligner), "", method).split()
                ]))

                curr_array = [round(len(fraction[fraction[method] == num]) *100 / len(fraction), 2)
                              for num in range(1, 6)]
                curr_array.append(round(
                    (len(fraction[fraction[method] == 6]) - found_in_all )*100 / len(fraction), 2))
                curr_array.append(round(found_in_all * 100 / len(fraction), 2))
                values_array.append(curr_array)

            plot.set_xticklabels(curr_labels, fontsize=10)
            values_array = numpy.array(values_array)
            values_array = values_array.transpose()
            X = numpy.arange(values_array.shape[1])
            for i in range(values_array.shape[0]):
                color = colors[(i+1)*2]
                bar = plot.bar(X, values_array[i],
                         bottom = numpy.sum(values_array[:i], axis=0),
                         color=color, edgecolor="grey")
                if row == col == 0:
                    # add handles to the legend
                    bar.get_children()[0].get_sketch_params()
                    legend_handles.append(bar)

    labels = ["Missed", "Intronic or Fragment", "Fusion", "Different structure",
                          "Contained", "Match", "Match in all"]
    plt.figlegend(labels=labels, framealpha=0.3,
                  loc="center right", handles=legend_handles,
                  fontsize=12,
                  ncol=1)
                  # ncol=floor(len(labels) / 2) + len(labels) % 2)

    plt.tight_layout(pad=0.15,
                     h_pad=1,
                     w_pad=0.5,
                     rect=[0.07,  # Left
                           0.1,  # Bottom
                           0.85,  # Right
                           0.9])  # Top

    if args.out is None:
        print("Showing the plot")
        plt.show()
    else:
        plt.savefig(args.out,
                    format=options["format"],
                    dpi=options["dpi"],
                    transparent=(not options["opaque"]))


def analyse_refmap(input_file, label, values):
    """
    Quick snippet to retrieve the ccodes from the RefMap
    :param input_file:
    :param label:
    :param values:
    :return:
    """

    with open(input_file) as refmap:
        ids_not_found = set()
        for row in csv.DictReader(refmap, delimiter="\t"):
            # "1.No Overlap",
            # "2. Fragment or intronic",
            # "3. Fusion",
            # "4. Alternative splicing",
            # "5. Extension (n, J)",
            # "6. Contained (c, C)",
            # "7. Match (=,_)"

            ccode = row["ccode"]
            if ccode in ("NA", "p", "P"):
                ccode = 1
            elif ccode in ("X", "x", "e", "I", "i", "ri", "rI"):
                ccode = 2
            elif ccode[0] == "f" and ccode.split(",")[1] not in ("=", "_"):
                ccode = 3
            elif ccode in ("G", "O", "g", "mo", "o", "h", "j"):
                ccode = 4
            elif ccode in ("n", "J"):
                ccode = 4
            elif ccode in ("c", "C"):
                ccode = 5
            else:
                assert (ccode in ("=", "_", "m") or (
                    ccode.split(",")[0] == "f" and ccode.split(",")[1] in ("=", "_"))), ccode
                ccode = 6

            if row["ref_id"] not in values:
                ids_not_found.add(row["ref_id"])
                continue
            values[row["ref_id"]][label] = ccode

        if len(ids_not_found) > 0:
            print("# of ids not found for {0}: {1}".format(
                input_file, len(ids_not_found)),
                  file=sys.stderr
            )

    return values


def main():

    parser = argparse.ArgumentParser("Script to create the merged ccode/TPM input for Gemy's modified script.")
    parser.add_argument("--quant_file", "-q",
                        required=True,
                        type=argparse.FileType("r"))
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--out", nargs="?", type=str, default=None,
                        help="Optional output file name. Default: None (show to stdout)")
    parser.add_argument("--title", type=str, default="")
    args = parser.parse_args()

    options = parse_configuration(args.configuration)
    values = dict()

    # Retrieve FPKM
    for row in csv.DictReader(args.quant_file, delimiter="\t"):
        values[row["target_id"]]=dict()
        values[row["target_id"]]["TPM"] = float(row["tpm"])
        values[row["target_id"]]["tid"] = row["target_id"]

    labels = []
    for aligner in sorted(options["divisions"]):
        for method in options["methods"]:
            label = "{} ({})".format(method, aligner)
            labels.append(label)
            input_file = "{}.refmap".format(
                re.sub(".stats", "", options["methods"][method][aligner][0]))
            values = analyse_refmap(input_file,
                                    label,
                                    values)

    tids = set(values.keys())
    right_total = len(labels) + 2
    count_zero = set()
    for tid in tids:
        if len(values[tid].keys()) == 2:
            del values[tid]
        elif values[tid]["TPM"] == 0:
            del values[tid]
            count_zero.add(tid)
        elif len(values[tid].keys()) != right_total:
            raise KeyError("ID {} has been found only in {}".format(tid, values[tid].keys()))
    tids = set.difference(tids, count_zero)
    if tids != set(values.keys()):
        print("Removed {} TIDs due to filtering".format(len(tids) - len(set(values.keys()))))

    data = defaultdict(list)

    sorter = functools.partial(sort_values, **{"dictionary": values})
    keys = None
    for tid in values:
        if keys is None:
            keys = values[tid].keys()
        for key in values[tid]:
            data[key].append(values[tid][key])

            #        out.writerow(values[tid])
    data = pandas.DataFrame(data, columns=["tid", "TPM"] + labels)
    generate_plot(data, args, options, nrows=len(options["divisions"]))

    return

main()
