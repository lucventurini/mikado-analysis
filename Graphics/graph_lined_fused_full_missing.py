import matplotlib.pyplot as plt
from collections import OrderedDict
from operator import itemgetter
import matplotlib.lines as mlines
import os
import csv
import argparse
from itertools import zip_longest
import yaml
import re


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("--log", action="store_true", default=False)
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    data = OrderedDict()

    options = yaml.load(args.configuration)
    args.configuration.close()

    for method in options["methods"]:
        for key in ("STAR", "colour", "index", "TopHat"):
            if key not in options["methods"][method]:
                raise KeyError("{} not found for {}".format(key.capitalize(),
                                                            method))
        for aligner in ("STAR", "TopHat"):
            if not isinstance(options["methods"][method][aligner], list):
                raise TypeError("Invalid type for aligner {}: {}".format(
                    aligner, type(options["methods"][method][aligner])))
            elif len(options["methods"][method][aligner]) != 2:
                raise ValueError("Invalid number of files specified for {} / {}".format(
                    method, aligner))
            elif any(not os.path.exists(_) for _ in options["methods"][method][aligner]):
                raise OSError("Files not found: {}".format(", ".join(
                    options["methods"][method][aligner])))

    new_methods = OrderedDict()

    for index, method in sorted([(options["methods"][method]["index"], method)
                          for method in options["methods"]], key=itemgetter(0)):
        new_methods[method] = options["methods"][method]

    options["methods"] = new_methods

    if len(set(options["methods"][method]["colour"]
               for method in options["methods"])) != len(options["methods"]):
        raise ValueError("Invalid unique number of colours specified!")

    if options["format"] not in ("svg", "png", "tiff"):
        raise ValueError("Invalid output format specified: {}".format(
            options["format"]))


    # text = """	Full	Missed	Fused
    # CLASS (STAR)	13428	7901	1702
    # CLASS (TopHat)	9493	10394	1300
    # Cufflinks (STAR)	11682	7494	2372
    # Cufflinks (TopHat)	10629	7172	2318
    # Stringtie (STAR)	13777	7154	2164
    # Stringtie (TopHat)	14755	6019	1985
    # Trinity (STAR)	14303	6113	870
    # Trinity (TopHat)	11266	7996	1310
    # Mikado (STAR)	16118	5988	633
    # Mikado (TopHat)	15721	6119	797
    # """

    for label in options["methods"]:
        for aligner in ("STAR", "TopHat"):
            data["{} ({})".format(label, aligner)] = [set(), set(), set()]
            orig_refmap = "{}.refmap".format(
                re.sub(".stats$", "", options["methods"][label][aligner][0]))
            with open(orig_refmap) as refmap:
                for row in csv.DictReader(refmap, delimiter="\t"):
                    if row["best_ccode"] in ("=", "_"):
                        data["{} ({})".format(label, aligner)][0].add(row["ref_gene"])
                    elif row["best_ccode"][0] == "f":
                        data["{} ({})".format(label, aligner)][2].add(row["ref_gene"])
                    elif row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                        data["{} ({})".format(label, aligner)][1].add(row["ref_gene"])
            for num in range(3):
                data["{} ({})".format(label, aligner)][num] = len(
                    data["{} ({})".format(label, aligner)][num])

    print(*data.items(), sep="\n")

    figure, axes = plt.subplots(nrows=1,
                                ncols=1,
                                dpi=300, figsize=(8, 6))
    # print(axes)
    plot = axes

    # Set the axis to log if necessary
    if args.log is True:
        plt.xscale("log")

    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (1, 1), 'k-')
    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (2, 2), 'k-')
    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (3, 3), 'k-')
    newticks = ["Fused genes", "Missed genes", "Recovered genes"]

    shape = ["o", "^"]

    handles = []
    for index, method in enumerate(data.keys()):
        color = options["methods"][method]["color"]
        # color = colors[int(index / 2)]
        marker = shape[index % 2]
        handle = mlines.Line2D([], [], markersize=5, color=color, marker=marker, label=method)
        # handles.append(handle)
        for pos, point in enumerate(reversed(data[method])):
            handle = plot.scatter(point, pos + 1, label=method, color=color, marker=marker, edgecolor="k", s=100)
            handle.get_sketch_params()
            if pos == 0:
                handles.append(handle)

    # handles, labels = plot.get_legend_handles_labels()
    plt.figlegend(labels=data.keys(), framealpha=0.3,
                  loc=(0.31, 0.09), handles=handles,
                  scatterpoints=1,
                  ncol=3, fontsize=10, markerscale=0.6)
    # plt.title("$Lines$")
    plot.set_ylim(0, 4)
    plot.set_xlim(0, max(max(data[_]) + 1000 for _ in data))
    plot.tick_params(axis='both', which='major', labelsize=10)
    plot.set_yticks([1, 2, 3])
    plot.set_yticklabels(newticks, fontsize=15)
    # plot.yticks([1, 2, 3], newticks, labelsize=8)
    plt.tight_layout(pad=0.5,
                     h_pad=1,
                     w_pad=1,
                     rect=[0.1,  # Left
                           0.25,  # Bottom
                           1,  # Right
                           1])  # Top
    if args.out is not None:
        plt.savefig(args.out, format="png", transparent=True)
    else:
        plt.show()

main()
