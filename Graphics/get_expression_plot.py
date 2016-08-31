#!/usr/bin/env python3

import sys
import argparse
import csv
import functools
import pandas
from collections import defaultdict
import matplotlib.cm as cm
import matplotlib.colors
# from rpy2.robjects import pandas2ri
import numpy
from matplotlib import pyplot as plt
# from matplotlib import ticker
# pandas2ri.activate()


def sort_values(key, dictionary):
    return dictionary[key]["TPM"]


def generate_plot(dataframe, args, nrows=2, ncols=4):

    # r = rpy2.robjects.r  # Start the R thread
    # base = importr("base")
    #
    # graphics = importr('graphics')
    #
    # # from rpy2.robjects.packages import importr
    # # from rpy2.robjects import pandas2ri
    # color_brewer = importr("RColorBrewer")
    # ggplot2 = importr("ggplot2")
    # grdevices = importr("grDevices")

    dataframe = dataframe.sort_values("TPM", ascending=True)
    print("Sorted dataframe", file=sys.stderr)
    greater_100 = dataframe[(dataframe.TPM > 10)]
    hundred_to_ten = dataframe[(dataframe.TPM <= 10) & (dataframe.TPM > 5)]
    ten_to_one = dataframe[(dataframe.TPM > 1) & (dataframe.TPM <= 5)]
    zer_to_one = dataframe[(dataframe.TPM > 0.01) & (dataframe.TPM <= 1)]
    lowest = dataframe[(dataframe.TPM <= 0.01)]
    print("Got the intervals", file=sys.stderr)

    color_map = cm.get_cmap(args.colourmap)
    color_normalizer = matplotlib.colors.Normalize(2, 8)

    colors = {1: "white"}

    for num in range(2, 8):
        colors[num] = color_map(color_normalizer(num))

    # plot = plt.plot()

    figure, axes=plt.subplots(nrows=nrows,
                              ncols=ncols,
                              figsize=(8, 6))

    methods = dataframe.columns[2:]
    current = 0

    newticks = ["<= 0.01", "0.01 - 1", "1-5", "5-10", "> 10"]

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
            _axes.set_xticks(numpy.arange(0.5,5.5,1))
            plot.set_xticklabels(newticks, fontsize=10)
            values_array = []
            for name, fraction in zip(
                    ("lowest", "zer_to_one", "ten_to_one", "hundred_to_ten", "greater_100"),
                    (lowest, zer_to_one, ten_to_one, hundred_to_ten, greater_100)):
                curr_array = [round(len(fraction[fraction[method] == num]) *100 / len(fraction), 2) for num in colors]
                values_array.append(curr_array)
            values_array = numpy.array(values_array)
            values_array = values_array.transpose()
            # print(method, values_array.shape, values_array)

            X = numpy.arange(values_array.shape[1])
            print(X)
            for i in range(values_array.shape[0]):
                plot.bar(X, values_array[i],
                         bottom = numpy.sum(values_array[:i], axis=0),
                         color=colors[i+1])
            plot.set_ylim(0, 100)
            plot.set_title("${}$".format(method))
            for tick in axes[row,col].get_xticklabels():
                tick.set_rotation(270)

    plt.tight_layout(pad=0.1,
                     h_pad=.2,
                     w_pad=0.2,
                     rect=[0.1,  # Left
                           0.2,  # Bottom
                           0.85,  # Right
                           0.9])  # Top

    # for curr_height, tid in enumerate(dataframe.tid[:100]):
    #     vals = dataframe[dataframe.tid == tid][dataframe.columns[2:]].astype(int).values[0]
    #     # print(curr_height, vals)
    #     for num in range(1,8):
    #         _bar = []
    #         for val in vals:
    #             if val == num:
    #                 _bar.append(1)
    #             else:
    #                 _bar.append(0)
    #         plt.bar(_bar, _bar, color=colors[num], bottom=curr_height)
    #     # print(curr_height)

    if args.out is None:
        plt.show()
    else:
        plt.savefig(args.out, format=args.format, dpi=args.dpi, transparent=args.opaque)
    # print("Creating the matrix")
    # data_matrix = base.data_matrix(dataframe[dataframe.columns[2:]])
    # print("Matrix created, setting the colnames")
    # data_matrix.colnames = rpy2.robjects.vectors.StrVector(dataframe.columns[2:])
    # print("Creating the output TIFF")
    #
    # tiff = grdevices.tiff(jpeg, width=13, height=13, units="in", res=300)
    # print("Created the output TIFF")
    #
    # colors = rpy2.robjects.vectors.StrVector(["white", "darkgrey", "black",
    #                                           # '"lightblue", "orange",
    #                                           "green", "red"])
    #
    # rowSide = rpy2.robjects.vectors.StrVector(["gray"]*greater_100 + ["violet"]*hundred_to_ten +
    #                                           ["sky blue"]*ten_to_one +
    #                                           ["light green"]*zer_to_one +
    #                                           ["yellow"]*lowest)
    # r.par(mai=rpy2.robjects.vectors.FloatVector([0, 0, 0, 0]))
    # print("Creating the HeatMap")
    # graphics.plot_new()
    # gplots = importr("gplots")
    # # print(base.dim(data_matrix))
    # # base.heatmap(data_matrix)
    # # gplots.heatmap_2(data_matrix, RowSideColors=rowSide)
    # # data_matrix = data_matrix[:1000]
    # gplots.heatmap_2(data_matrix,
    #                  RowSideColors=rowSide,
    #                  col=colors,
    #                  dendogram="none",
    #                  breaks=numpy.arange(0.5, 8.5),
    #                  density_info="none",
    #                  trace="none",
    #                  margins=rpy2.robjects.vectors.IntVector([18,15]),
    #                  cexCol=2, cexMain=10,
    #                  # par=r.par(cex_main=10),
    #                  cex_axis=1, cex_main=5, cex_lab=3, cex_sub=3,
    #                  col_axis="red", col_lab="red",
    #                  key=False,
    #                  Rowv="NA",
    #                  Colv="NA")
    #
    # #                  density_info="none",
    #
    # print("Creating the legend")
    # r.legend("top", fill=colors,
    #          legend=rpy2.robjects.vectors.StrVector(["1.No Overlap",
    #                                                  "2. Fragment or intronic",
    #                                                  "3. Fusion",
    #                                                  "4. Alternative splicing",
    #                                                  "5. Match (=,_)"]),
    #          cex=1.8,
    #          title="Class code legend")
    #
    #
    # r.legend("left",
    #          legend=rpy2.robjects.vectors.StrVector([">100 TPM", "10-100 TPM",
    #                                                  "1-10 TPM", "0.01-1 TPM", "upto 0.01 TPM"]),
    #          col=rpy2.robjects.vectors.StrVector(["gray", "violet",
    #                                               "sky blue", "light green", "yellow"]),
    #          cex=1.8,
    #          lty=1, lwd=10, title="TPM")
    #
    # print("Finished")
    #
    # grdevices.dev_off()


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
                ccode = 5
            elif ccode in ("c", "C"):
                ccode = 6
            else:
                assert (ccode in ("=", "_", "m") or (
                    ccode.split(",")[0] == "f" and ccode.split(",")[1] in ("=", "_"))), ccode
                ccode = 7

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
    parser.add_argument("-l", "--labels", default=None, type=str, nargs="+",
                        help="Labels for the input files, comma separated. Required.", required=True)
    parser.add_argument("-cm", "--colourmap", default="Accent",
                        help="Colourmap to be used.")
    parser.add_argument("--dpi", default=600, type=int)
    parser.add_argument("--format", default="svg", choices=["svg",
                                                            "pdf",
                                                            "png"])
    parser.add_argument("--out", nargs="?", type=str, default=None,
                        help="Optional output file name. Default: None (show to stdout)")
    parser.add_argument("input_files", help="The RefMap input files", nargs="+")
    args = parser.parse_args()

    # args.labels = args.labels.split(",")
    if len(args.labels) != len(args.input_files):
        raise ValueError("Labels must be the same number as input files!")        

    # Create data dictionary
    values = dict()

    # Retrieve FPKM
    for row in csv.DictReader(args.quant_file, delimiter="\t"):
        values[row["target_id"]]=dict()
        values[row["target_id"]]["TPM"] = float(row["tpm"])
        values[row["target_id"]]["tid"] = row["target_id"]

    for input_file, label in zip(args.input_files, args.labels):
        values = analyse_refmap(input_file, label, values)

    data = defaultdict(list)

    sorter = functools.partial(sort_values, **{"dictionary": values})
    keys = None
    for tid in values:
        if keys is None:
            print(values[tid].keys())
            keys = values[tid].keys()
        for key in values[tid]:
            data[key].append(values[tid][key])

            #        out.writerow(values[tid])
    data = pandas.DataFrame(data, columns=["tid", "TPM"] + args.labels)
    generate_plot(data, args)

    return

main()
