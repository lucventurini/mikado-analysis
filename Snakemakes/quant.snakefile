import os
import glob


rule all:
    input: os.path.join("Reference", "reference.cdnas.kallisto"), os.path.join("Reference", "reference.cdnas.fasta.fai")
    output: touch("quant.done")

rule kallisto_index:
    input: os.path.join("Reference", "reference.cdnas.fasta")
    output: os.path.join("Reference", "reference.cdnas.kallisto")
    message: "set +u && source kallisto-0.43.0 && kallisto index -i {output} {input} && set -u"
    shell: "set +u && source kallisto-0.43.0 && kallisto index -i {output} {input} && set -u"

rule faidx:
    input:
        os.path.join("Reference", "reference.cdnas.fasta")
    output:
        os.path.join("Reference", "reference.cdnas.fasta.fai")
    shell: "set +u && source samtools-1.2 && samtools faidx {input}"

rule extract:
    input:
        ref=os.path.join("Reference", "reference.gff3"),
	genome=os.path.join("Reference", "genome.fa")
    output:
        os.path.join("Reference", "reference.cdnas.fasta")
    shell: "set + u && source cufflinks-2.1.1 && gffread -w {output} -g {input.genome} {input.ref}"