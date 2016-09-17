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


def sort_values(key, dictionary):
    return dictionary[key]["TPM"]


def generate_plot(dataframe, args, options, nrows=2, ncols=5):

    dataframe = dataframe.sort_values("TPM", ascending=True)
    greater_100 = dataframe[(dataframe.TPM > 100)]
    hundred_to_ten = dataframe[(dataframe.TPM <= 100) & (dataframe.TPM > 10)]
    ten_to_five = dataframe[(dataframe.TPM > 5) & (dataframe.TPM <= 10)]
    five_to_one = dataframe[(dataframe.TPM > 1) & (dataframe.TPM <= 5)]
    zer_to_one = dataframe[(dataframe.TPM > 0.01) & (dataframe.TPM <= 1)]
    lowest = dataframe[(dataframe.TPM <= 0.01)]
    plt.style.context("ggplot")

    color_map = cm.get_cmap(options["colourmap"]["name"])
    color_normalizer = matplotlib.colors.Normalize(2, 15)

    colors = {1: "white", 2: "white"}

    for num in range(3, 7*2 + 1):
        colors[num] = color_map(color_normalizer(num))

    # plot = plt.plot()

    figure, axes=plt.subplots(nrows=nrows,
                              ncols=ncols,
                              figsize=(16, 10))

    methods = dataframe.columns[2:]
    current = 0

    newticks = ["0 - 0.01", "0.01 - 1", "1-5", "5-10", "10-100", "$\geq$ 100"]

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
    figure.text(0.02, 0.6, "$Transcripts$ $per$ $category$ $(\%)$", va="center", fontsize=15, rotation="vertical")

    legend_handles = []
    for row in range(nrows):
        for col in range(ncols):
            if current == len(methods):
                axes[row, col].xaxis.set_visible(False)
                axes[row, col].yaxis.set_visible(False)
                continue
            method = methods[current]
            current +=1
            plot = axes[row, col]
            _axes = plot.axes
            # _axes.xaxis.set_major_locator(ticker.MultipleLocator(1))
            _axes.set_xticks(numpy.arange(0.5, 6.5,1))
            plot.set_xticklabels(newticks, fontsize=10)
            values_array = []
            for name, fraction in zip(
                    ("lowest", "zer_to_one", "five_to_one", "ten_to_five", "hundred_to_ten", "greater_100"),
                    (lowest, zer_to_one, five_to_one, ten_to_five, hundred_to_ten, greater_100)):
                curr_array = [round(len(fraction[fraction[method] == num]) *100 / len(fraction), 2) for num in range(1, 7)]
                values_array.append(curr_array)

            values_array = numpy.array(values_array)
            values_array = values_array.transpose()
            # print(method, values_array.shape, values_array)

            X = numpy.arange(values_array.shape[1])
            for i in range(values_array.shape[0]):
                bar = plot.bar(X, values_array[i],
                         bottom = numpy.sum(values_array[:i], axis=0),
                         color=colors[(i+1)*2])
                if row == col == 0:
                    # add handles to the legend
                    bar.get_children()[0].get_sketch_params()
                    legend_handles.append(bar)

            plot.set_ylim(0, 100)
            plot.set_title(" ".join(["${}$".format(_) for _ in method.split()]))
            for tick in axes[row,col].get_xticklabels():
                tick.set_rotation(60)

    labels = ["Missed", "Intronic or Fragment", "Fusion", "Different structure",
                          "Contained", "Match"]
    plt.figlegend(labels=labels, framealpha=0.3,
                  loc="lower center", handles=legend_handles,
                  ncol=floor(len(labels) / 2) + len(labels) % 2)

    plt.tight_layout(pad=0.15,
                     h_pad=.2,
                     w_pad=0.1,
                     rect=[0.1,  # Left
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

    options = parse_configuration(args)
    values = dict()

    # Retrieve FPKM
    for row in csv.DictReader(args.quant_file, delimiter="\t"):
        values[row["target_id"]]=dict()
        values[row["target_id"]]["TPM"] = float(row["tpm"])
        values[row["target_id"]]["tid"] = row["target_id"]

    labels = []
    for aligner in ("STAR", "TopHat"):
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
    for tid in tids:
        if len(values[tid].keys()) == 2:
            del values[tid]
        elif values[tid]["TPM"] == 0:
            del values[tid]
        elif len(values[tid].keys()) != right_total:
            raise KeyError("ID {} has been found only in {}".format(tid, values[tid].keys()))
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
    generate_plot(data, args, options)

    return

main()
