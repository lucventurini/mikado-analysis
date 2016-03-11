#!/usr/bin/env python3

from Mikado.parsers.GFF import GFF3
import csv
import argparse
import functools

__author__ = 'Luca Venturini'


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("loci_refmap", type=argparse.FileType("rt"))
    parser.add_argument("subloci_refmap", type=argparse.FileType("rt"))
    parser.add_argument("loci_gff", type=GFF3)
    args = parser.parse_args()

    loci_refmap = csv.DictReader(args.loci_refmap, delimiter="\t")
    sub_refmap = csv.DictReader(args.subloci_refmap, delimiter="\t")

    loci_res = dict()

    alias_dict = dict()

    for row in args.loci_gff:
        if row.is_transcript is False:
            continue
        alias_dict[row.id] = row.attributes["alias"]

    args.loci_gff.close()

    for row in loci_refmap:
        if row["tid"] != "NA":
            row["tid"] = alias_dict[row["tid"]]
        loci_res[row["ref_id"]] = row

    for row in sub_refmap:
        if row["ccode"] in ("=", "_") and loci_res[row["ref_id"]]["ccode"] not in ("=","_"):
            print("Failed:", row["tid"], row["ccode"], loci_res[row["ref_id"]]["tid"], loci_res[row["ref_id"]]["ccode"])
        elif row["ccode"] not in ("=", "_") and loci_res[row["ref_id"]]["ccode"] in ("=","_"):
            print("Recovered:", row["tid"], row["ccode"], loci_res[row["ref_id"]]["tid"], loci_res[row["ref_id"]]["ccode"])

    args.loci_refmap.close()
    args.subloci_refmap.close()

main()

