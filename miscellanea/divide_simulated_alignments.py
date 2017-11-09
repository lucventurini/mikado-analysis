#!/urs/bin/env python3


import pysam
import argparse
import sys
import re
import os


__doc__ = """"""

def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("bam", type=pysam.AlignmentFile)
    args = parser.parse_args()

    mismapped = pysam.AlignmentFile(args.out + ".mismapped.bam", "wb",
                                    header=args.bam.header)
    correct = pysam.AlignmentFile(args.out + ".correct.bam", "wb",
                                  header=args.bam.header)

    pattern = re.compile(".*_([^:]*):([0-9]*):([^_]*)_([0-9]*)_([^#]*).*")

    for record in args.bam:

        if record.mate_is_unmapped is True or record.is_unmapped is True:
            mismapped.write(record)
            continue

        chrom, pos, next_chrom, next_pos, template = (record.reference_name, record.reference_start,
                                                      record.next_reference_name, record.next_reference_start,
                                                      record.tlen)

        pat_chrom, pat_pos, pat_fcigar, pat_template, pat_scigar = re.search(pattern, record.qname).groups()
        pat_pos = int(pat_pos) - 1
        # pat_template = int(pat_template)

        if chrom != pat_chrom or next_chrom != pat_chrom:
            mismapped.write(record)
            continue
        elif ((record.is_read1 and pat_fcigar != record.cigarstring)
              or (record.is_read2 and pat_scigar != record.cigarstring)):
            mismapped.write(record)
            continue
        elif pat_pos not in (record.reference_start, record.next_reference_start):
            mismapped.write(record)
            continue
        else:
            correct.write(record)

    mismapped.close()
    correct.close()

main()