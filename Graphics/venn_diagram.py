# import readline
# import rpy2.robjects
import csv
import multiprocessing
from collections import OrderedDict, Counter
# from rpy2.robjects.packages import importr
import itertools
import argparse
import matplotlib.cm as cm
import matplotlib.colors
import intervene.modules.venn.list_venn as ivenn
import re
import sys
from utils import parse_configuration, parse_refmaps
# import numpy


def clamp(x):
    return max(0, min(x, 255))

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
                        required=True, nargs="+")
    parser.add_argument("--transcripts", action="store_true", default=False,
                        help="Flag. If set, Venn plotted against transcripts, not genes.")
    parser.add_argument("--dpi", default=300, type=int)
    parser.add_argument("--title", default="Venn Diagram")
    parser.add_argument("--procs", default=1, type=int)
    args = parser.parse_args()

    options = parse_configuration(args.configuration, exclude_mikado=args.exclude)

    sets = OrderedDict.fromkeys(["{}\n({})".format(*_) for _ in itertools.product(options["methods"], args.aligner)])

    for k in list(sets.keys()) + ["base"]:
        sets[k] = {"full": set(), "fusion": set(), "missing": set()}

    # Update the sets for each gene and label
    if args.transcripts is True:
        colname = "ref_id"
        ccode = "ccode"
        tag = "transcripts"
    else:
        colname = "ref_gene"
        ccode = "best_ccode"
        tag = "genes"

    first = True

    pool = multiprocessing.Pool(processes=args.procs)
    proxies = dict().fromkeys(sets.keys())

    for method in options["methods"]:
        for aligner in args.aligner:
            if options["methods"][method][aligner] is not None:
                orig_stats, filtered_stats = options["methods"][method][aligner][:2]
            else:
                orig_stats, filtered_stats = None, None
            if first is True and orig_stats is not None:
                proxies["base"] = pool.apply_async(parse_refmaps, (orig_stats, filtered_stats, args.transcripts, True, True))
            key = "{}\n({})".format(method, aligner)
            proxies[key] = pool.apply_async(parse_refmaps, (orig_stats, filtered_stats, args.transcripts, True, False))

    for proxy in proxies:
        sets[proxy]["full"], sets[proxy]["missing"], sets[proxy]["fusion"] = proxies[proxy].get()

    print("Loaded RefMaps.", file=sys.stderr)

    # Now use intervene venn

    labels = dict()
    sums = dict()
    for typ in ["full", "missing", "fusion"]:
        labels[typ] = ivenn.get_labels([list(sets[_][typ]) for _ in sets])
        sums[typ] = dict()
        for num in range(len(sets) + 1):
            sums[typ][num] = 0
        for label in labels[typ]:
            # print(typ, label, sum(int(_) for _ in label), int(labels[typ][label]))
            sums[typ][sum(int(_) for _ in label) - 1] += int(labels[typ][label])
            continue

    print("Sums")
    for num in sorted(range(len(sets))):
        print(num, *[sums[_][num] for _ in ["full", "missing", "fusion"]])

    print("Per method")
    all_reconstructable = set.union(*[sets[_]["full"] for _ in sets if _ != "base"])
    missed_all = set.intersection(*[sets[_]["missing"] for _ in sets if _ != "base"])
    fused_all = set.intersection(*[sets[_]["fusion"] for _ in sets if _ != "base"])

    for ds in sets:
        if ds == "base":
            continue
        key = " ".join(ds.split("\n"))
        row = [key]
        row.append(len(sets[ds]["full"]))
        row.append(len(set.difference(all_reconstructable, sets[ds]["full"])))
        row.append(round(100 * len(set.difference(all_reconstructable, sets[ds]["full"]))/len(all_reconstructable), 2))
        print(*row, sep="\t")
        
    # print("Labels:", labels[args.type])

    if len(sets) > 6:
        print("Too many sets to intersect ({}), exiting.".format(len(sets)), file=sys.stderr)
        sys.exit(0)

    funcs = {2: ivenn.venn2,
             3: ivenn.venn3,
             4: ivenn.venn4,
             5: ivenn.venn5,
             6: ivenn.venn6,}

    # Recalculate labels without the base
    labels = dict()
    for typ in ["full", "missing", "fusion"]:
        labels[typ] = ivenn.get_labels([list(sets[_][typ]) for _ in sets if _ != "base"])

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
                nums = (int(matched.groups()[0]), int(matched.groups()[1]), int(matched.groups()[2]))
                if nums == (255, 255, 255):  # Pure white
                    nums = (125, 125, 125)
                cols[index] = "#{0:02x}{1:02x}{2:02x}{3:02x}".format(clamp(nums[0]), clamp(nums[1]), clamp(nums[2]), 80)
            
    print(labels[args.type])
    print([_ for _ in sets.keys() if _ != "base"])
    print(cols)
    if (len(sets) - 1) == 2:
        fontsize = 18
    else:
        fontsize = 20
    fig, ax = funcs[len(sets) - 1](labels[args.type],
                               names=[_ for _ in sets.keys() if _ != "base"],
                               colors=cols,
                               fontsize=fontsize,
                               dpi=args.dpi,
                               alpha=0.5,
                               figsize=(12, 12))
    fig.savefig("{}.{}".format(args.out, args.format),
                dpi=args.dpi)
    print("Saved the figure to {}.{}".format(args.out, args.format))
    import time
    time.sleep(3)

main()
