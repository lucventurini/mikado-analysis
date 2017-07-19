#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
import argparse
import Mikado
from collections import defaultdict
import io

def main():

    parser = argparse.ArgumentParser("""Script to define the reference genes reconstructable using the input RNA-Seq data, for the Mikado publication.""")
    parser.add_argument("-r", "--reference", type=argparse.FileType("r"), required=True,
                        help="Reference annotation in GFF3 format")
    parser.add_argument("-p", "--portcullis", type=argparse.FileType("r"), required=True,
                        help="Junctions present in the RNA-Seq data, derived using Portcullis, in GFF3 format")
    parser.add_argument("-b", "--bam", type=argparse.FileType("r"), required=True,
                        help="Input BAM file.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output prefix")
    parser.add_argument("-s", "--strand-specific", dest="strand_specific", action="store_true",
                        default=False)
    args = parser.parse_args()

    args.bam.close()
    args.bam = args.bam.name

    args.portcullis.close()
    args.portcullis = args.portcullis.name

    args.reference.close()
    args.reference = args.reference.name    
    
    logger = Mikado.utilities.log_utils.create_default_logger("inference")
    logger.setLevel("INFO")

    logger.info("""Starting to analyse the dataset.
Commandline: {cli}
BAM: {bam}
Reference: {ref}
Junctions: {juncs}""".format(cli=" ".join(sys.argv),
                             bam=args.bam,
                             ref=args.reference,
                             juncs=args.portcullis))
    
    # Check the presence of required program
    gt = shutil.which("gt")
    if gt is None:
        logger.critical("Genome tools not found! Please install or source before continuing")
        sys.exit(1)
    else:
        logger.info("Gt binary: {0}".format(gt))

    bedtools = shutil.which("bedtools")
    if bedtools is None:
        logger.critical("Bedtools not found! Please install or source before continuing")
        sys.exit(1)
    else:
        logger.info("Bedtools binary: {0}".format(bedtools))

    # Derive the exon GFF
    logger.info("Starting to print out the exon GFF3")
    count = 0
    tids = dict()

    exons = []
    if not os.path.exists("exons.gff"):
        with open("exons.gff", "w") as exon_out, open(args.reference) as ref:
            for line in Mikado.parsers.GFF.GFF3(ref):
                if line.feature == "exon":
                    count += 1
                    exons.append(line)
                elif line.is_transcript is True:
                    tids[line.id] = line.parent
        print(*sorted(exons), sep="\n", file=exon_out)
        pass
    else:
        with open(args.reference) as ref:
            for line in Mikado.parsers.GFF.GFF3(ref):
                if line.feature == "exon":
                    count += 1
                    exons.append(line)
                elif line.is_transcript is True:
                    tids[line.id] = line.parent
                    
                    
    if count == 0:
        logger.critical("No exon found in the reference annotation!")
        sys.exit(1)
    else:
        logger.info("Found {0} exons for {1} transcripts".format(count, len(tids)))

    # Derive the introns
    logger.info("Starting to create the intron GFF")
    intron_creator = subprocess.Popen(["gt", "gff3", "-addintrons", "-retainids", args.reference],
                                      shell=False, stdout=subprocess.PIPE)
    count = 0
    reference_introns = defaultdict(list)
    multiexonic_tids = set()
    for line in intron_creator.stdout:
        line = line.decode()
        line = Mikado.parsers.GFF.GffLine(line)
        if line.feature == "intron":
            key = (line.chrom, line.start, line.end)  # We do not keep the strand because it might be a non-canonical one
            reference_introns[key].append(line.parent[0])  # Store the name of the parent
            multiexonic_tids.add(line.parent[0])
            count += 1

    # Get all the TIDs of transcripts which do not have introns
            
    if count == 0:
        logger.critical("No intron found in the reference! Exiting")
        sys.exit(1)
    elif  intron_creator.returncode not in (None, 0):
        logger.critical("Genometools malfunctioned, exiting")
        sys.exit(1)
    else:
        logger.info("Found {0} introns in the reference".format(count))

    monoexonic_tids = set.difference(set(tids.keys()), multiexonic_tids)
    logger.info("Out of {0} transcripts, {1} are monoexonic ({2}%) and {3} are multiexonic ({4}%)".format(
        len(tids), len(monoexonic_tids), round(100 * len(monoexonic_tids) / len(tids), 2),
        len(multiexonic_tids), round(100 * len(multiexonic_tids) / len(tids), 2)))

    found_introns = set()
    with open(args.portcullis) as portcullis:
        for line in Mikado.parsers.GFF.GFF3(portcullis):
            key = (line.chrom, line.start, line.end)
            found_introns.add(key)

    tids_not_intron_covered = set()
    for key in reference_introns:
        if key not in found_introns:
            tids_not_intron_covered.update(reference_introns[key])
    
    logger.info("{0} transcripts had at least one intron not covered by this RNA-Seq data.".format(len(tids_not_intron_covered)))
    tids_all_intron_covered = set.difference(multiexonic_tids, tids_not_intron_covered)
    logger.info("{0} multiexonic transcripts had all their introns covered".format(len(tids_all_intron_covered)))

    # Use Bedtools to derive the coverageBed
    bedtools_cli = ["bedtools", "coverage", "-split", "-sorted", "-a", "exons.gff", "-b", args.bam]
    if args.strand_specific is True:
        bedtools_cli.append("-s")

    if not os.path.exists("coverage_exons.gff"):
        logger.info("Deriving the bedtools coverage for the exons, commandline:\n\t\t\t{0}".format(" ".join(bedtools_cli)))

        with open("coverage_exons.gff", "w") as coverage_exon_out:
            bedtools_command = subprocess.call(bedtools_cli, shell=False, stdout=coverage_exon_out)
        
        if bedtools_command not in (0, None):
            logger.critical("Something has gone wrong with the coverage exon command!")
            sys.exit(1)
        else:
            logger.info("Finished the bedtools coverage command")
    else:
        logger.info("coverage_exons.gff already created, skipping BEDtools")
        with open("coverage_exons.gff") as coverage_exon_out:
            # This is just to set the coverage_exons_out handle variable.
            pass

    with_uncovered_exons = set()
    with open(coverage_exon_out.name) as coverage_exon_out:
        for line in coverage_exon_out:
            fields = line.rstrip().split("\t")
            gffline, coverage = Mikado.parsers.GFF.GffLine("\t".join(fields[:9])), float(fields[-1])
            if coverage < 1:
                with_uncovered_exons.add(gffline.parent[0])

    logger.info("{0} transcripts had at least one exon not completely covered".format(len(with_uncovered_exons)))
    reconstructable_multi = set.difference(tids_all_intron_covered, with_uncovered_exons)
    reconstructable_mono = set.difference(monoexonic_tids, with_uncovered_exons)
    reconstructable_all = set.union(reconstructable_multi, reconstructable_mono)
    if len(reconstructable_all) != len(reconstructable_mono) + len(reconstructable_multi):
        logger.critical("Something has gone awry, total number of recs ({0}) is different from sum of recs multi and mono ({1} and {2})".format(
            len(reconstructable_all), len(reconstructable_mono), len(reconstructable_multi)))
        sys.exit(1)
    
    logger.info("{0} multiexonic and {1} monoexonic transcripts reconstructable (total {2})".format(
        len(reconstructable_multi), len(reconstructable_mono), len(reconstructable_all)))

    with open(args.output,'w') as out:
        for tid in sorted(list(reconstructable_all)):
            print(tid, tids[tid], file=out, sep="\t")

if __name__ == "__main__":
    main()
                
    
                        
