import rpy2.robjects
import sys
import csv
import os
import rpy2.robjects
from collections import OrderedDict, Counter
from rpy2.robjects.packages import importr
import itertools
import argparse
import matplotlib.cm as cm
import matplotlib.colors


def main():
    
    parser = argparse.ArgumentParser("Script to create the Venn Plots")
    parser.add_argument("-t",
                        "--type",
                        choices=["missing", "full", "fusion"],
                        required=True)
    parser.add_argument("-r", "--refmap", required=True,
                        nargs="+")
    parser.add_argument("-l", "--labels", required=True,
                        nargs="+")
    parser.add_argument("-cm", "--colourmap", default="gist_rainbow",
                        help="Colourmap to be used.")
    parser.add_argument("-c", "--colours", nargs="+",
                        help="Colors to use.")
    parser.add_argument("-o", "--out",
                        type=str, help="Output file", default="venn.svg")
    parser.add_argument("--format", choices=["svg", "tiff", "png"], default="tiff")
    parser.add_argument("--transcripts", action="store_true", default=False,
                        help="Flag. If set, Venn plotted against transcripts, not genes.")
    parser.add_argument("--title", default="Venn Diagram")
    args = parser.parse_args()

    if len(args.refmap) != len(args.labels):
        print("Incorrect number of labels specified!")
        parser.print_help()
        sys.exit(1)
    elif len(args.refmap) > 5:
        print("It is impossible for me to create a Venn plot with more than 5 circles.")
        parser.print_help()
        sys.exit(1)

    sets = OrderedDict.fromkeys(args.labels)
    for k in sets:
        sets[k] = set()

    total = Counter()
    first = True

    # Update the sets for each gene and label
    if args.transcripts is True:
        colname = "ref_id"
        tag="transcripts"
    else:
        colname = "ref_gene"
        tag = "genes"
    for val, refmap in zip(args.labels, args.refmap):
        tsv = csv.DictReader(open("{0}".format(refmap)), delimiter="\t")
        for row in tsv:
            if first:
                total.update([row[colname]])
            if row["ccode"].lower() in ("na", "x", "p", "i", "ri") and args.type == "missing":
                sets[val].add(row[colname])
            elif row["ccode"] in ("=", "_") and args.type == "full":
                sets[val].add(row[colname])
            elif row["ccode"][0] == "f" and args.type == "fusion":
                sets[val].add(row[colname])
            else:
                continue
        if first:
            for gid in total:
                total[gid] = 0
            first = False

    r = rpy2.robjects.r  # Start the R thread                                                                                       
    base = importr("base")
    venn = importr("VennDiagram")
    grdevices = importr("grDevices")
    corrs = dict((x+1, list(sets.keys())[x]) for x in range(len(sets.keys())))
    nums = dict()

    for num in corrs:
        cat = corrs[num]
        nums["area{0}".format(num)] = len(sets[cat])
        total.update(list(sets[cat]))
        print(cat.capitalize(), nums["area{0}".format(num)])

    print("Total", len(set.union(*sets.values())))
    counts = list(total.values())

    for num in range(len(args.labels), 0, -1):
        tot = counts.count(num)
        cum_tot = sum(counts.count(x) for x in range(num+1, 6)) + tot
        print("Genes {0} by {1} methods ({2} cumulative): {3}".format(
            args.type, num, cum_tot, tot))
    print("")

    if len(args.labels) == 2:
        nums["cross.area"] = len(set.intersection(
            sets[corrs[1]], sets[corrs[2]]))
    else:
        for num_combs in range(2, len(args.labels) + 1):
            for comb in itertools.combinations(range(1, len(args.labels) + 1), num_combs):
                # print(comb)
                index = "".join([str(x) for x in comb])
                curr_sets = [sets[corrs[num]] for num in comb]
                nums["n{0}".format(index)] = len(set.intersection(*curr_sets))

    # for num_combs in range(2,6):
    #     for comb in itertools.combinations(range(1,6), num_combs):
    #         index = "".join([str(x) for x in comb])
    #         curr_sets = [sets[corrs[num]] for num in comb]
    #         nums["n{0}".format(index)] = len(set.intersection(*curr_sets))

    # Get the colours from the ColorMap
    # return
    if args.colours is None:
        color_normalizer = matplotlib.colors.Normalize(0, len(args.labels))
        color_map = cm.get_cmap(args.colourmap)
        cols = [matplotlib.colors.rgb2hex(color_map(color_normalizer(index))) for index in range(len(args.labels))]
        cols = rpy2.robjects.vectors.StrVector(cols)
    else:
        print(args.colours)
        cols = rpy2.robjects.vectors.StrVector(args.colours)
    #

    if len(args.labels) == 1:
        draw_function = venn.draw_single_venn
    elif len(args.labels) == 2:
        draw_function = venn.draw_pairwise_venn
    elif len(args.labels) == 3:
        draw_function = venn.draw_triple_venn
    elif len(args.labels) == 4:
        draw_function = venn.draw_quad_venn
    else:
        draw_function = venn.draw_quintuple_venn

    if not args.transcripts:
        distances = {
            1: [0.1],
            2: [0.1, 0.1],
            3: [0.1, 0.1, 0.1],
            4: [0.3, 0.25, 0.15, 0.15],
            5: [0.2, 0.3, 0.2, 0.25, 0.25]
        }
    else:
        distances = {
            1: [0.1],
            2: [0.1, 0.1],
            3: [0.1, 0.1, 0.1],
            4: [0.35, 0.3, 0.2, 0.25],
            5: [0.2, 0.3, 0.2, 0.25, 0.25]
        }

    # draw_function = venn.venn_diagram
    gridExtra = importr("gridExtra")
    grid = importr("grid")
    dev_args = {"width": 960, "height": 960}
    if args.format == "tiff":
        device = grdevices.tiff
    elif args.format == "png":
        device = grdevices.png
        dev_args["bg"] = "transparent"
    else:
        device = grdevices.svg
    
    device(args.out, **dev_args)

    drawn = draw_function(height=4000, width=4000,
                          fill=cols,
                          category=rpy2.robjects.vectors.StrVector(["{}\n({} {})".format(x.capitalize(), len(sets[x]), tag) for x in sets.keys()]),
                          margin=0.2,
                          cat_dist=rpy2.robjects.vectors.FloatVector(distances[len(args.labels)]),
                          cat_cex=3,
                          cat_col=cols,
                          cex=2,
                          cex_main=args.title,
                          # main_pos=rpy2.robjects.vectors.FloatVector([2500, 100]),
                          **nums)
    gridExtra.grid_arrange(grid.gTree(children=drawn),
                           top=grid.textGrob(args.title, gp=grid.gpar(fontsize=30)))
    grdevices.dev_off()

main()
