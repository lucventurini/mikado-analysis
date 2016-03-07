#!/usr/bin/env python3

import sys
from matplotlib import pyplot as pl
import argparse


def main():

    # Command line:
    # /usr/local/bin/Mikado.py compare -r reference.gff3 -p Mikado.py.loci.gff3 \
    # -o compare -l compare.log
    # 7 reference RNAs in 5 genes
    # 22 predicted RNAs in  14 genes
    # ------------------------------ |   Sn |   Pr |   F1 |
    #                      Base level: 85.79  50.20  63.33
    #          Exon level (stringent): 63.83  36.59  46.51
    #            Exon level (lenient): 80.00  49.32  61.02
    #                    Intron level: 89.47  58.62  70.83
    #              Intron chain level: 33.33  12.50  18.18
    #    Transcript level (stringent): 0.00  0.00  0.00
    #      Transcript level (lenient): 28.57  9.09  13.79
    #          Gene level (stringent): 0.00  0.00  0.00
    #            Gene level (lenient): 40.00  14.29  21.05

    #          Matching intron chains: 2
    #           Matched intron chains: 2

    #        Missed exons (stringent): 17/47  (36.17%)
    #         Novel exons (stringent): 52/82  (63.41%)
    #          Missed exons (lenient): 9/45  (20.00%)
    #           Novel exons (lenient): 37/73  (50.68%)
    #                  Missed introns: 4/38  (10.53%)
    #                   Novel introns: 24/58  (41.38%)

    #              Missed transcripts: 0/7  (0.00%)
    #               Novel transcripts: 9/22  (40.91%)
    #                    Missed genes: 0/5  (0.00%)
    #                     Novel genes: 5/14  (35.71%)

    parser = argparse.ArgumentParser("Script to generate the plots for the Mikado.py stats")
    parser.add_argument("stat", type=argparse.FileType("r"),
                        nargs="+",
                        help="The Mikado stats file")
    args = parser.parse_args()

    levels = {"base": dict(),
              "exon": dict(),
              "intron": dict(),
              "intron chain": dict(),
              "transcript": dict(),
              "gene": dict()}

    for stat_file in args.stat:
        for line in stat_file:
            fields = line.strip().split()
            level = None
            if len(fields) < 3:
                continue
            if fields[0] == "Base":
                level = "base"
            elif fields[0] == "Intron":
                if fields[1] == "chain":
                    level = "intron chain"
                else:
                    level = "intron"
            elif fields[2] == "(lenient):":
                if fields[0] == "Exon":
                    level = "exon"
                elif fields[0] == "Transcript":
                    level = "transcript"
                elif fields[0] == "Gene":
                    level = "gene"
                else:
                    print(fields[0])
                    continue
            else:
                print([fields[0], fields[2]])
            if level is None:
                continue
            f1 = float(fields[-1])
            levels[level][stat_file.name] = f1

            pass
    print(levels)

    return

if __name__ == '__main__':
    main()