import sys
# from utils import line_correspondence, parse_configuration
from utils import line_correspondence, parse_configuration
import argparse
import tabulate
import yaml
import numpy as np
from math import ceil
from itertools import zip_longest
from collections import OrderedDict
from scipy.stats.mstats import hmean
import warnings


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def calc_f1(rcl, prc):

    if rcl == 0 or prc == 0:
        return 0
    else:
        return 2*(rcl*prc)/(rcl+prc)


def rcl(f1s, prc):
    denominator = (2.0 * prc - f1s)
    if denominator == 0:
        return 100
    else:
        return f1s * prc / denominator


def main():

    """"""

    parser = argparse.ArgumentParser(__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--species", type=argparse.FileType("r"), required=True)
    parser.add_argument("--out", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("--format", choices=tabulate._table_formats.keys(), default="latex")
    parser.add_argument("--level", default="transcript", choices=line_correspondence.keys())
    args = parser.parse_args()

    species = yaml.load(args.species)

    # Names are the SPECIES
    names = species.pop("names")
    # Categories are Real vs Simulated data
    categories = species.pop("categories")

    # We want to create a table of the form

    # Species | Method   || Category              || Category              || ..
    # species | Name     | Prec | Rec | F1        ||Prec | Rec | F1        || ..

    name_ar = []
    for category in categories:
        for name in names:
            # print(category, name)
            key = [_ for _ in species.keys() if
                   isinstance(species[_], dict) and species[_]["name"] == name and species[_]["category"] == category]
            # print(key)
            assert len(key) == 1
            key = key.pop()
            name_ar.append(key)
    name_ar = np.array(list(grouper(name_ar, ceil(len(species) / 2), None)))

    header = [""] + list(categories)

    header.append(["Species", "Aligner", "Method"] + ["Precision", "Recall", "F1"] * 2)

    rows = []
    # print(rows)

    key = None
    methods = None
    divisions = None

    for yrow, name in enumerate(names):

        new_rows = OrderedDict()

        for xrow, category in enumerate(categories):
            try:
                key = name_ar[xrow, yrow]
            except IndexError:
                raise IndexError(name_ar, xrow, yrow)
            if key is None:
                raise IndexError(name_ar, xrow, category, yrow, name)
                # continue

            with open(species[key]["configuration"]) as configuration:
                options = parse_configuration(configuration, prefix=species[key]["folder"])
                # Assembler
                if methods is None:
                    methods = list(options["methods"])
                    divisions = list(options["divisions"])

                for method in options["methods"]:
                    # Aligner
                    for division in options["divisions"]:
                        meth_key = (method, division)
                        if meth_key not in new_rows:
                            new_rows[meth_key] = OrderedDict()
                        try:
                            orig, filtered = options["methods"][method][division]
                        except TypeError:
                            warnings.warn("Something went wrong for {}, {}; continuing".format(
                                method, division))

                            new_rows[meth_key][category] = (-10, -10, -10)
                            continue
                        orig_lines = [line.rstrip() for line in open(orig)]
                        filtered_lines = [line.rstrip() for line in open(filtered)]
                        for index, line_index in enumerate([line_correspondence[args.level]]):
                            precision = float(orig_lines[line_index].split(":")[1].split()[1])
                            recall = float(filtered_lines[line_index].split(":")[1].split()[0])
                            try:
                                f1 = hmean(np.array([precision, recall]))
                            except TypeError as exc:
                                raise TypeError("\n".join([str(_) for _ in [(precision, type(precision)),
                                                                            (recall, type(recall)),
                                                                            exc]]))
                            # print(level, method, division, (precision, recall, f1))
                            new_rows[meth_key][category] = (precision, recall, f1)

        begun = False
        for division in divisions:
            division_done = False
            for method in methods:
                meth_key = (method, division)
                if not begun:
                    row = [name]
                    begun = True
                else:
                    row = [""]

                if not division_done:
                    row.append(division)
                    division_done = True
                else:
                    row.append("")

                row.append(method)
                # row.append(meth_key)
                # print(new_rows[meth_key].keys())
                for category in new_rows[meth_key]:
                    row.extend(new_rows[meth_key][category])

                rows.append(row)

    print(tabulate.tabulate(rows, headers=header, tablefmt=args.format))
    # print(categories)
    return



if __name__ == "__main__":
    main()