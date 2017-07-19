import readline
import rpy2.robjects
import csv
from collections import OrderedDict, Counter
from rpy2.robjects.packages import importr
import itertools
import argparse
import matplotlib.cm as cm
import matplotlib.colors
import intervene.modules.venn.list_venn as ivenn
import re
from utils import parse_configuration
import numpy


def clamp(x):
    return max(0, min(x, 255))


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
    parser.add_argument("--format", choices=["svg", "tiff", "png"], default="svg")
    parser.add_argument("-a", "--aligner", # choices=["STAR", "TopHat"],
                        required=True)
    parser.add_argument("--transcripts", action="store_true", default=False,
                        help="Flag. If set, Venn plotted against transcripts, not genes.")
    parser.add_argument("--dpi", default=300, type=int)
    parser.add_argument("--title", default="Venn Diagram")
    args = parser.parse_args()

    options = parse_configuration(args.configuration, exclude_mikado=args.exclude)

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
    print("Loaded RefMaps.")

    # Now use intervene venn

    labels = ivenn.get_labels([list(sets[_]) for _ in sets])
    print("Labels:", labels)
    funcs = {2: ivenn.venn2,
             3: ivenn.venn3,
             4: ivenn.venn4,
             5: ivenn.venn5,
             6: ivenn.venn6,}

    if options["colourmap"]["use"] is True:
        color_normalizer = matplotlib.colors.Normalize(0, len(options["methods"]))
        color_map = cm.get_cmap(options["colourmap"]["name"])
        cols = [color_map(color_normalizer(index))
                for index in range(len(options["methods"]))]
        # cols = [matplotlib.colors.rgb2hex(color_map(color_normalizer(index)))
        #         for index in range(len(options["methods"]))]
        # cols = rpy2.robjects.vectors.StrVector(cols)
    else:
        cols = [options["methods"][_]["colour"] for _ in options["methods"]]
        for index, colour in enumerate(cols):
            matched = re.match("\(([0-9]*), ([0-9]*), ([0-9]*)\)$", colour)
            if matched:
                colour = "#{0:02x}{1:02x}{2:02x}{3:02x}".format(clamp(int(matched.groups()[0])),
                                                                 clamp(int(matched.groups()[1])),
                                                                 clamp(int(matched.groups()[2])),
                                                                50)

            cols[index] = colour
            
        print(cols)
        # cols = rpy2.robjects.vectors.StrVector(cols)

    fig, ax = funcs[len(sets)](labels, names=list(options["methods"].keys()),
                               colors=cols,
                               fontsize=100,
                               dpi=args.dpi,
                               alpha=0.5,
                               figsize=(7, 7) if len(options) < 5 else (12, 12))
    fig.savefig("{}.{}".format(args.out, args.format),
                dpi=args.dpi)
    print("Saved the figure to {}.{}".format(args.out, args.format))
    import time
    time.sleep(3)

main()
