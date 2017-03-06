#!/usr/bin/env python3

import sys
import argparse
import csv
# import functools
import re
import pandas
# from collections import defaultdict
# import matplotlib.cm as cm
# import matplotlib.colors
import numpy
# import scipy.stats
import matplotlib.pyplot as plt
# import matplotlib.patches
# from math import floor
import seaborn
from utils import parse_configuration
import os


def sort_values(key, dictionary):
    return dictionary[key]["TPM"]


# def reject_outliers(data):
#     extreme = scipy.stats.mstats.mquantiles(data, prob=[.01, .99])
#     data[]
#
#     d = numpy.abs(data - numpy.median(data))
#     mdev = numpy.median(d)
#     s = d/mdev if mdev else 0.
#     return data[s<m]


def generate_plot(star_data, tophat_data, merged_data, typ, args, options):

    plt.style.context("ggplot")
    print("Sorting values")
    star_data = star_data.sort_values(["TPM", "tid"], ascending=True)
    tophat_data = tophat_data.sort_values(["TPM", "tid"], ascending=True)
    merged_data = merged_data.sort_values(["TPM", "tid"], ascending=True)
    # Remove outliers
    print("Sorted values")

    print("Prepared the data")

    figure, axes = plt.subplots(nrows=3,
                                ncols=1,
                                figsize=(16, 10),
                                sharex=False,
                                sharey=True)

    figure.suptitle(" ".join(["${}$".format(re.sub("%", "\%", _)) for _ in args.title.split()]),
                    fontsize=20, style="italic", family="serif")

    # Plot complete genes
    if typ == "full":
        x_num = 6
        line = "Number of methods capable of correctly reconstructing the transcript"
    elif typ == "missed":
        x_num = 1
        line = "Number of methods which have completely missed the transcript"
    else:
        raise ValueError("Invalid type: {}".format(typ))

    for index, (method, dataframe, color) in enumerate([("STAR", star_data, "turquoise"),
                                                        ("TopHat", tophat_data, "coral"),
                                                        ("Either", merged_data, "slategrey")]):
        plt.axes(axes[index])
        seaborn.violinplot(data=dataframe,
                           color=color,
                           # hue="Aligner",
                           # palette="Set2",
                           x=x_num,
                           y="TPM",
                           cut=2,
                           inner="box",
                           showmeans=True,
                           scale="width")
        x_max = max(dataframe[x_num].max() + 1, dataframe[x_num].max() + 1)
        axes[index].set_xlim(-1, x_max)
        axes[index].set_xlabel(line + " ({})".format(method), size=18)
        axes[index].set_ylabel("$log_{2}(TPM+1)$", size=20)
        # axes[index].xaxis.set_major_locator(plt.NullLocator())
        # for tick in axes[0].xaxis.get_major_ticks():
        #     tick.label.set_fontsize(0)
        for tick in axes[index].yaxis.get_major_ticks():
            tick.label.set_fontsize(13)
        for tick in axes[index].yaxis.get_major_ticks():
                tick.label.set_fontsize(13)

    plt.tight_layout(pad=0.2,
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
        out = os.path.splitext(args.out)[0] + ".{}".format(typ) + os.path.splitext(args.out)[1] + ".".format(args.format)
        plt.savefig(out,
                    format=options["format"],
                    dpi=options["dpi"],
                    transparent=(not options["opaque"]))


def analyse_refmap(input_file, values):
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
            values[row["ref_id"]][ccode] += 1

        if len(ids_not_found) > 0:
            print("# of ids not found for {0}: {1}".format(
                input_file, len(ids_not_found)),
                  file=sys.stderr
            )

    return values


def main():

    parser = argparse.ArgumentParser("Script to create a violin plot.")
    parser.add_argument("--quant_file", "-q",
                        required=True,
                        type=argparse.FileType("r"))
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--out", nargs="?", type=str, default=None,
                        help="Optional output file name. Default: None (show to stdout)")
    parser.add_argument("--format", choices=["svg", "png", "pdf"],
                        default="svg")
    parser.add_argument("--title", type=str, default="")
    args = parser.parse_args()

    options = parse_configuration(args)
    # Retrieve FPKM
    star = dict()
    tophat = dict()

    for row in csv.DictReader(args.quant_file, delimiter="\t"):
        # Set up STAR
        star[row["target_id"]]=dict()
        star[row["target_id"]]["TPM"] = numpy.log2(float(row["tpm"])+1)
        star[row["target_id"]]["tid"] = row["target_id"]
        star[row["target_id"]]["Aligner"] = "STAR"
        for num in range(1, 7):
            star[row["target_id"]][num] = 0

        # Set up TopHat
        tophat[row["target_id"]] = dict()
        tophat[row["target_id"]]["TPM"] = numpy.log2(float(row["tpm"])+1)
        tophat[row["target_id"]]["tid"] = row["target_id"]
        tophat[row["target_id"]]["Aligner"] = "TopHat"
        for num in range(1, 7):
            tophat[row["target_id"]][num] = 0

    for aligner, dictionary in (("STAR", star), ("TopHat", tophat)):
        for method in options["methods"]:
            input_file = "{}.refmap".format(
                re.sub(".stats", "", options["methods"][method][aligner][0]))
            print(input_file)
            dictionary = analyse_refmap(input_file,
                                        dictionary)

    star_data = []
    tophat_data = []
    merged_data = []
    for tid in star:
        star_data.append([tid] + [star[tid]["TPM"]] + ["STAR"] +
                    [star[tid][_] for _ in range(1,7)])
        tophat_data.append([tid] + [tophat[tid]["TPM"]] + ["Tophat"] +
                    [tophat[tid][_] for _ in range(1, 7)])
        merged_data.append([tid] + [tophat[tid]["TPM"]] + ["Either"] +
                           [star[tid][_] + tophat[tid][_] for _ in range(1, 7)])

    print("Generating the final data frames")
    star_data = pandas.DataFrame(star_data, columns=["tid", "TPM", "Aligner"] + list(range(1,7)))
    tophat_data = pandas.DataFrame(tophat_data, columns=["tid", "TPM", "Aligner"] + list(range(1,7)))
    merged_data = pandas.DataFrame(merged_data, columns=["tid", "TPM", "Aligner"] + list(range(1,7)))

    generate_plot(star_data, tophat_data, merged_data, "full", args, options)
    generate_plot(star_data, tophat_data, merged_data, "missed", args, options)
    return

main()
