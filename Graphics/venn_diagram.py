import rpy2.robjects
import csv
import rpy2.robjects
from collections import OrderedDict, Counter
from rpy2.robjects.packages import importr
import itertools
import argparse
import matplotlib.cm as cm
import matplotlib.colors
import re
from .utils import parse_configuration


def main():
    
    parser = argparse.ArgumentParser("Script to create the Venn Plots")
    parser.add_argument("-t",
                        "--type",
                        choices=["missing", "full", "fusion"],
                        required=True)
    parser.add_argument("-c", "--configuration", required=True, type=argparse.FileType("r"))
    parser.add_argument("-em", "--exclude-mikado", dest="exclude",
                        action="store_true", default=False,
                        help="Flag. If set, Mikado results will be excluded")
    parser.add_argument("-o", "--out",
                        type=str, help="Output file", required=True)
    parser.add_argument("--format", choices=["svg", "tiff", "png"], default=None)
    parser.add_argument("-a", "--aligner", choices=["STAR", "TopHat"],
                        required=True)
    parser.add_argument("--transcripts", action="store_true", default=False,
                        help="Flag. If set, Venn plotted against transcripts, not genes.")
    parser.add_argument("--title", default="Venn Diagram")
    args = parser.parse_args()

    options = parse_configuration(args, exclude_mikado=args.exclude)

    sets = OrderedDict.fromkeys(options["methods"])

    for k in sets:
        sets[k] = set()

    total = Counter()
    first = True

    # Update the sets for each gene and label
    if args.transcripts is True:
        colname = "ref_id"
        ccode = "ccode"
        tag = "transcripts"
    else:
        colname = "ref_gene"
        ccode = "best_ccode"
        tag = "genes"

    for method in options["methods"]:
        refmap = "{}.refmap".format(
            re.sub(".stats$", "", options["methods"][method][args.aligner][0]))
        with open(refmap) as ref:
            tsv = csv.DictReader(ref, delimiter="\t")
            for row in tsv:
                if first:
                    total.update([row[colname]])
                if row[ccode].lower() in ("na", "x", "p", "i", "ri") and args.type == "missing":
                    sets[method].add(row[colname])
                elif row[ccode] in ("=", "_") and args.type == "full":
                    sets[method].add(row[colname])
                elif row[ccode][0] == "f" and args.type == "fusion":
                    sets[method].add(row[colname])
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

    for num in range(len(options["methods"]), 0, -1):
        tot = counts.count(num)
        cum_tot = sum(counts.count(x) for x in range(num+1, 6)) + tot
        print("Genes {0} by {1} methods ({2} cumulative): {3}".format(
            args.type, num, cum_tot, tot))
    print("")

    if len(options["methods"]) == 2:
        nums["cross.area"] = len(set.intersection(
            sets[corrs[1]], sets[corrs[2]]))
    else:
        for num_combs in range(2, len(options["methods"]) + 1):
            for comb in itertools.combinations(range(1, len(options["methods"]) + 1), num_combs):
                # print(comb)
                index = "".join([str(x) for x in comb])
                curr_sets = [sets[corrs[num]] for num in comb]
                nums["n{0}".format(index)] = len(set.intersection(*curr_sets))

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(args.colourmap)
        cols = [matplotlib.colors.rgb2hex(color_map(color_normalizer(index)))
                for index in range(len(options["methods"]))]
        cols = rpy2.robjects.vectors.StrVector(cols)
    else:
        cols = [options["methods"][_]["colour"] for _ in options["methods"]]
        cols = rpy2.robjects.vectors.StrVector(cols)

    if len(options["methods"]) == 1:
        draw_function = venn.draw_single_venn
    elif len(options["methods"]) == 2:
        draw_function = venn.draw_pairwise_venn
    elif len(options["methods"]) == 3:
        draw_function = venn.draw_triple_venn
    elif len(options["methods"]) == 4:
        draw_function = venn.draw_quad_venn
    else:
        draw_function = venn.draw_quintuple_venn

    if not args.transcripts:
        distances = {
            1: [0.1],
            2: [0.1, 0.1],
            3: [0.1, 0.1, 0.1],
            4: [0.3, 0.25, 0.15, 0.15],
            5: [0.25, 0.37, 0.25, 0.27, 0.33]
        }
    else:
        distances = {
            1: [0.1],
            2: [0.1, 0.1],
            3: [0.1, 0.1, 0.1],
            4: [0.35, 0.3, 0.2, 0.25],
            5: [0.3, 0.3, 0.27, 0.28, 0.28]
        }

    # draw_function = venn.venn_diagram
    gridExtra = importr("gridExtra")
    grid = importr("grid")
    dev_args = {"width": 960, "height": 960}
    if options["format"] == "tiff":
        device = grdevices.tiff
    elif options["format"] == "png":
        device = grdevices.png
        dev_args["bg"] = "transparent"
    else:
        device = grdevices.svg
    
    device(args.out, **dev_args)

    drawn = draw_function(height=4000, width=4000,
                          fill=cols,
                          category=rpy2.robjects.vectors.StrVector(
                              ["{}\n({:,}\n {})".format(
                                  "\n".join([_.capitalize() for _ in x.split()]),
                                  len(sets[x]), tag) for x in sets.keys()]),
                          margin=0.2,
                          cat_dist=rpy2.robjects.vectors.FloatVector(distances[
                                                                         len(options["methods"])]),
                          cat_cex=2.5,
                          cat_col=cols,
                          cex=2,
                          cex_main=args.title,
                          # main_pos=rpy2.robjects.vectors.FloatVector([2500, 100]),
                          **nums)
    gridExtra.grid_arrange(grid.gTree(children=drawn),
                           top=grid.textGrob(args.title, gp=grid.gpar(fontsize=30)))
    grdevices.dev_off()

main()
