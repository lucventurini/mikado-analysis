#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import functools
import pandas
from collections import defaultdict
import rpy2.robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import numpy

pandas2ri.activate()


def sort_values(key, dictionary):
    return dictionary[key]["TPM"]


def generate_plot(dataframe, jpeg):

    r = rpy2.robjects.r  # Start the R thread
    base = importr("base")
    
    graphics = importr('graphics')

    # from rpy2.robjects.packages import importr
    # from rpy2.robjects import pandas2ri
    color_brewer = importr("RColorBrewer")
    ggplot2 = importr("ggplot2")
    grdevices = importr("grDevices")

    dataframe = dataframe.sort_values("TPM", ascending=False)
    print("Sorted dataframe", file=sys.stderr)
    greater_100 = len(dataframe[(dataframe.TPM > 100)].tid)
    hundred_to_ten = len(dataframe[(dataframe.TPM <= 100) & (dataframe.TPM > 10)].tid)
    ten_to_one = len(dataframe[(dataframe.TPM > 1) & (dataframe.TPM <= 10)].tid)
    zer_to_one = len(dataframe[(dataframe.TPM > 0.01) & (dataframe.TPM <= 1)].tid)
    lowest = len(dataframe[(dataframe.TPM <= 0.01)].tid)
    print("Got the intervals", file=sys.stderr)

    print("Creating the matrix")
    data_matrix = base.data_matrix(dataframe[dataframe.columns[2:]])
    print("Matrix created, setting the colnames")
    data_matrix.colnames = rpy2.robjects.vectors.StrVector(dataframe.columns[2:])
    print("Creating the output TIFF")

    tiff = grdevices.tiff(jpeg, width=13, height=13, units="in", res=300)
    print("Created the output TIFF")
    
    colors = rpy2.robjects.vectors.StrVector(["white", "darkgrey", "black",
                                              "lightblue", "orange", "green", "red"])

    rowSide = rpy2.robjects.vectors.StrVector(["gray"]*greater_100 + ["violet"]*hundred_to_ten +
                                              ["sky blue"]*ten_to_one +
                                              ["light green"]*zer_to_one +
                                              ["yellow"]*lowest)
    r.par(mai=rpy2.robjects.vectors.FloatVector([0, 0, 0, 0]))
    print("Creating the HeatMap")
    graphics.plot_new()
    gplots = importr("gplots")
    gplots.heatmap_2(data_matrix, RowSideColors=rowSide)
    # gplots.heatmap_2(data_matrix, RowSideColors=rowSide,
    #                  density_info="none",
    #                  trace="none",
    #                  margins=rpy2.robjects.vectors.IntVector([18,15]),
    #                  col=colors,
    #                  breaks=numpy.arange(0.5, 6.5),
    #                  cexCol=2, cexMain=10,
    #                  # par=r.par(cex_main=10),
    #                  cex_axis=1, cex_main=5, cex_lab=3, cex_sub=3,
    #                  col_axis="red", col_lab="red",
    #                  dendogram="none", key=False,
    #                  Rowv="NA", Colv="NA")
    print("Creating the legend")

    r.legend("top", fill=colors,
             legend=rpy2.robjects.vectors.StrVector(["1.No Overlap",
                                                     "2. Fragment or intronic",
                                                     "3. Fusion",
                                                     "4. Alternative splicing",
                                                     "5. Extension (n, J)",
                                                     "6. Contained (c, C, m)",
                                                     "7. Match (=,_)"]),
             cex=1.8,
             title="Class code legend")

    r.legend("left",
             legend = rpy2.robjects.vectors.StrVector([">100 TPM", "10-100 TPM",
                                                       "1-10 TPM", "0.01-1 TPM", "upto 0.01 TPM"]),
             col = rpy2.robjects.vectors.StrVector(["gray", "violet",
                                                    "sky blue", "light green", "yellow"]),
             cex=1.8,
             lty=1, lwd=10, title="TPM")
    print("Finished")

    grdevices.dev_off()


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
    parser.add_argument("-l", "--labels", default=None, type=str,
                        help="Labels for the input files, comma separated. Required.", required=True)
    parser.add_argument("--out", nargs="?", type=str, default="expression.tiff",
                        help="Optional output file name. Default: %(default)s")
    parser.add_argument("input_files", help="The RefMap input files", nargs="+")
    args = parser.parse_args()

    args.labels = args.labels.split(",")
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
    generate_plot(data, args.out)

    return

main()
