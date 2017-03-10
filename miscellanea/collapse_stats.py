#!/usr/bin/env python3

import argparse
import sys
import re

__doc__ = "Script to collapse various stat files into one."

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", default=sys.stdout,
                        type=argparse.FileType("w"))
    parser.add_argument("stat", type=argparse.FileType("rt"),
                        nargs="+")
    args = parser.parse_args()

    data = []
    # 1 Command line:
    # 2 /tgac/software/testing/bin/core/../..//mikado/1.0.0b1/x86_64/bin/mikado compare -r Reference/reference.gff3 -p Assemblies/STAR/Trinity/results.gtf -l Comparisons/STAR/Original/Trinity.log -o Comparisons/STAR/Original/Trinity
    # 3 35386 reference RNAs in 27416 genes
    # 4 105212 predicted RNAs in  44798 genes
    # 17
    # 18 #   Matching: in prediction; matched: in reference.
    # 19
    # 26
    # 33
    
    for stat in args.stat:
        new = dict()
        new["stat"] = []
        new["matches"] = []
        new["features"] = []
        new["high_level"] = []
        for num, line in enumerate(stat, start=1):
            line = line.rstrip()
            # 1 Command line:
            # 2 /tgac/software/testing/mikado/0.19.0/x86_64/bin/mikado.py compare [..]
            # 3 27643 reference RNAs in 14475 genes
            # 4 16142 predicted RNAs in  11379 genes
            if num < 3:
                continue
            elif num == 3:
                ref_transcr, ref_genes = [int(_) for _ in re.sub(" reference RNAs in", "", re.sub(" genes", "", line)).split()]
                new["ref_transcr"] = ["{0},,".format(ref_transcr)]
                new["ref_genes"] = ["{0},,".format(ref_genes)]
            elif num == 4:
                pred_transcr, pred_genes = [int(_) for _ in re.sub(" predicted RNAs in", "",
                                                                 re.sub(" genes", "", line)).split()]
                new["pred_transcr"] = ["{0},,".format(pred_transcr)]
                new["pred_genes"] = ["{0},,".format(pred_genes)]
            elif num in range(6, 17):
                # 5 --------------------------------- |   Sn |   Pr |   F1 |
                # 6                         Base level: 74.90  56.67  64.52
                # 7             Exon level (stringent): 49.97  34.23  40.63
                # 8               Exon level (lenient): 73.45  65.81  69.42
                # 9                       Intron level: 78.48  87.81  82.89
                # 10                 Intron chain level: 45.87  23.68  31.24
                # 11       Transcript level (stringent): 0.01  0.00  0.01
                # 12   Transcript level (>=95% base F1): 25.94  10.29  14.73
                # 13   Transcript level (>=80% base F1): 41.03  18.32  25.33
                # 14          Gene level (100% base F1): 0.01  0.01  0.01
                # 15         Gene level (>=95% base F1): 31.82  19.63  24.28
                # 16         Gene level (>=80% base F1): 49.80  31.96  38.94
                new["stat"].append(",".join(line.split(":")[1].lstrip().split()))
            elif num in range(20, 26):
                # 20             Matching intron chains: 18473
                # 21              Matched intron chains: 13354
                # 22    Matching monoexonic transcripts: 2098
                # 23     Matched monoexonic transcripts: 2096
                # 24         Total matching transcripts: 20571
                # 25          Total matched transcripts: 15450
                new["matches"].append("{0},,".format(line.split(":")[1].lstrip()))
            elif num in range(27, 33):
                # 27           Missed exons (stringent): 79866/159647  (50.03%)
                # 28            Novel exons (stringent): 153284/233065  (65.77%)
                # 29             Missed exons (lenient): 40068/150903  (26.55%)
                # 30              Novel exons (lenient): 57581/168416  (34.19%)
                # 31                     Missed introns: 26803/124567  (21.52%)
                # 32                      Novel introns: 13568/111332  (12.19%)
                new["features"].append("{0},,".format(line.split(":")[1].lstrip()))
            elif num in range(34,38):
                # 34                 Missed transcripts: 5658/35386  (15.99%)
                # 35                  Novel transcripts: 1116/105212  (1.06%)
                # 36                       Missed genes: 5316/27416  (19.39%)
                # 37                        Novel genes: 761/44798  (1.70%)
                new["high_level"].append("{0},,".format(line.split(":")[1].lstrip()))

            else:
                continue
        data.append(new)
        continue

    for key in ["stat", "matches", "features", "high_level",
                "ref_transcr", "ref_genes", "pred_transcr", "pred_genes"]:
        for row in zip(*[_[key] for _ in data]):
            print(",".join(row), file=args.out)
            continue
        if not ("pred" in key or "ref" in key):
            print("", file=args.out)
        else:
            pass
        continue
    
    # assert len(data) == len(args.stat), len(data)
    # print(data)
            
main()
