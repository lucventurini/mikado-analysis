#!/usr/bin/env python3

import argparse
import sys
from collections import Counter, defaultdict
import pandas

__doc__ = """This quick script will output, given a set of refmaps, how many transcripts
have been correctly reconstructed by all of them, all minus one, all minus two, etc."""


codes = {'=': 6,
         'C': 5,
         'G': 4,
         'I': 2,
         'J': 4,
         'NA': 1,
         'O': 4,
         'P': 1,
         'X': 2,
         '_': 6,
         'c': 5,
         'e': 2,
         'f': 3,
         'g': 4,
         'h': 4,
         'i': 2,
         'j': 4,
         'm': 6,
         'mo': 4,
         'n': 4,
         'o': 4,
         'p': 1,
         'rI': 2,
         'ri': 2,
         'x': 2}


def analyse_refmap(input_file, values, feature="transcript"):
    """
    Quick snippet to retrieve the ccodes from the RefMap
    :param input_file:
    :param label:
    :param values:
    :return:
    """

    if feature == "transcript":
        ccode_key = "ccode"
        unique = False
        key = "ref_id"
    else:
        ccode_key = "best_ccode"
        unique = True
        key = "ref_gene"

    df = pandas.read_csv(input_file, delimiter="\t")
    cols = df[[key, ccode_key]]
    if unique is True:
        cols = cols.drop_duplicates()

    cols.fillna("NA")
    cols[ccode_key].replace("f,=", "=", inplace=True, regex=True)
    cols[ccode_key].replace("f,_", "_", inplace=True, regex=True)
    cols[ccode_key].replace("f,.*", "f", inplace=True, regex=True)
    cols[ccode_key] = cols[ccode_key].map(codes)

    for row in cols.itertuples():
        _, rid, val = row
        if rid not in values:
            values[rid] = defaultdict(int)
        values[rid][val] += 1

    return values


def main():

    keys = {"full": 6,
            "missed": 1,
            "fragment": 2,
            "fusion": 3,
            "as": 4,
            "contained": 5}

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("-t", "--type", default="full",
                        choices=keys.keys())
    parser.add_argument("-f", "--feature", choices=["transcript", "gene"], default="gene")
    parser.add_argument("refmap", nargs="+")
    args = parser.parse_args()

    vals = dict()
    [analyse_refmap(refmap, vals, feature=args.feature) for refmap in args.refmap]

    counter = Counter(vals.keys())
    key = keys[args.type]
    for ref_id in vals:
        counter.update([vals[ref_id][key]])

    for num in reversed(range(0, len(args.refmap) + 1)):
        print(num, counter[num], round(100 * counter[num] / len(vals), 2), sep="\t", file=args.out)

main()



