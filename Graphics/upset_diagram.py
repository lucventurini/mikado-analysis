# import rpy2.robjects
import csv
# import rpy2.robjects
from collections import OrderedDict, Counter
# from rpy2.robjects.packages import importr
# import itertools
import argparse
# import matplotlib.cm as cm
# import matplotlib.colors
import re
from utils import parse_configuration
import pandas as pd
import pyupset as pyu
import matplotlib.pyplot as plt


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
    # parser.add_argument("-a", "--aligner", choices=["STAR", "TopHat"],
    #                     required=True)
    parser.add_argument("--transcripts", action="store_true", default=False,
                        help="Flag. If set, Venn plotted against transcripts, not genes.")
    parser.add_argument("--title", default="Venn Diagram")
    args = parser.parse_args()

    options = parse_configuration(args, exclude_mikado=args.exclude)

    sets = OrderedDict()

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

    for aligner in ["STAR", "TopHat"]:
        for method in options["methods"]:
            refmap = "{}.refmap".format(
                re.sub(".stats$", "", options["methods"][method][aligner][0]))
            with open(refmap) as ref:
                tsv = csv.DictReader(ref, delimiter="\t")
                meth = "{} ({})".format(method, aligner)
                sets[meth] = set()
                for row in tsv:
                    if first:
                        total.update([row[colname]])
                    if row[ccode].lower() in ("na", "x", "p", "i", "ri") and args.type == "missing":
                        sets[meth].add(row[colname])
                    elif row[ccode] in ("=", "_") and args.type == "full":
                        sets[meth].add(row[colname])
                    elif row[ccode][0] == "f" and args.type == "fusion":
                        sets[meth].add(row[colname])
                    else:
                        continue
                if first:
                    for gid in total:
                        total[gid] = 0
                    first = False

    for aligner in ["STAR", "TopHat"]:
        for method in sorted(options["methods"].keys()):
            set_name = "{} ({})".format(method, aligner)
            print(set_name)
            sets[set_name] = pd.DataFrame(list(sets[set_name]), columns=["TID"])

    pyu.plot(sets,
             # sort_by="degree",
             inters_size_bounds=(100, 20000),
             )
    if args.format is None:
        args.format = "svg"
    plt.savefig(args.out, format=args.format)

main()
