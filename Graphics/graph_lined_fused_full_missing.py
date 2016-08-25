import matplotlib.pyplot as plt
from collections import OrderedDict as odict
import matplotlib.lines as mlines
import os
import csv
import argparse
from itertools import zip_longest


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--labels", nargs=5, required=True)
    parser.add_argument("-o", "--out", default=None)
    parser.add_argument("--star", nargs=5)
    parser.add_argument("--tophat", nargs=5)
    # parser.add_argument("refmap", nargs=10, type=argparse.FileType("rt"))
    args = parser.parse_args()

    data = odict()

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

    for label, star, th in zip(args.labels, args.star, args.tophat):
        for aligner, orig_refmap in zip(["STAR", "TopHat"], [star, th]):
            data["{} ({})".format(label, aligner)] = [set(), set(), set()]
            with open(orig_refmap) as refmap:
                for row in csv.DictReader(refmap, delimiter="\t"):
                    if row["best_ccode"] in ("=", "_"):
                        data["{} ({})".format(label, aligner)][0].add(row["ref_gene"])
                    elif row["best_ccode"][0] == "f":
                        data["{} ({})".format(label, aligner)][2].add(row["ref_gene"])
                    elif row["best_ccode"] in ("NA", "p", "P", "i", "I", "ri", "rI", "X", "x"):
                        data["{} ({})".format(label, aligner)][1].add(row["ref_gene"])
            # with open(filter_refmap) as refmap:
            #     for row in csv.DictReader(refmap, delimiter="\t"):

            for num in range(3):
                data["{} ({})".format(label, aligner)][num] = len(
                    data["{} ({})".format(label, aligner)][num])

    print(*data.items(), sep="\n")

    # for row in text.split("\n"):
    #     if row.startswith("\t"):
    #         continue
    #     aligner, *points = row.split("\t")
    #     if len(points) == 0:
    #         continue
    #     data[aligner] = [int(_) for _ in points]
    #     continue

    figure, axes = plt.subplots(nrows=1,
                                ncols=1,
                                dpi=300, figsize=(8, 6))
    # print(axes)
    plot = axes

    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (1, 1), 'k-')
    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (2, 2), 'k-')
    plot.plot((1, max(max(data[_]) + 1000 for _ in data)), (3, 3), 'k-')
    newticks = ["Fused genes", "Missed genes", "Recovered genes"]

    colors = ["lightgreen", "pink", "blue", "yellow", "red"]
    shape = ["o", "^"]

    handles = []
    for index, method in enumerate(data.keys()):
        color = colors[int(index / 2)]
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
