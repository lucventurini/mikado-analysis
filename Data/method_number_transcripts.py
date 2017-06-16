#!/usr/bin/env python3

import csv
import argparse
import sys
from collections import defaultdict

__doc__ = """This quick script will output, given a set of refmaps, how many transcripts
have been correctly reconstrcted by all of them, all minus one, all minus two, etc."""


def analyse_refmap(input_file, values):
    """
    Quick snippet to retrieve the ccodes from the RefMap
    :param input_file:
    :param values:
    :return:
    """

    with open(input_file) as refmap:
        # ids_not_found = set()
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
                ccode = 4
            elif ccode in ("c", "C"):
                ccode = 5
            else:
                assert (ccode in ("=", "_", "m") or (
                    ccode.split(",")[0] == "f" and ccode.split(",")[1] in ("=", "_"))), ccode
                ccode = 6
            #
            # if row["ref_id"] not in values:
            #     values[row["ref_id"]] = defaultdict(int)

            # if row["ref_id"] not in values:
            #     ids_not_found.add(row["ref_id"])
            #     continue
            values[ccode] += 1

        # if len(ids_not_found) > 0:
        #     print("# of ids not found for {0}: {1}".format(
        #         input_file, len(ids_not_found)),
        #           file=sys.stderr
        #     )

    return values


def main():

    keys = {"full": 6,
            "missed": 1,
            "fragment": 2,
            "fusion": 3,
            "as": 4,
            "extension": 5}

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("-t", "--type", default="full",
                        choices=keys.keys())
    parser.add_argument("txt")
    args = parser.parse_args()

    with open(args.txt) as txt:
        for line in txt:
            tag, fname = line.rstrip().split("\t")
            vals = analyse_refmap(fname, defaultdict(int))
            print(tag, vals[keys[args.type]],
                  round(vals[keys[args.type]]/sum(vals.values()), 4),
                  sep="\t", file=args.out
                  )

main()
