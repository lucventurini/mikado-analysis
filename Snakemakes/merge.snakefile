import sys
import glob
import subprocess
import os

aln_methods = ["STAR", "TopHat"]
asm_methods = ["CLASS", "Cufflinks", "Stringtie", "Trinity"]

rule all:
     # input: expand(os.path.join("Assemblies", "{aligner}", "{assembler}", "results.gtf"), aligner=aln_methods, assembler=asm_methods)
     # output:
     input:
         merged = expand(os.path.join("CuffMerge", "{aligner}", "Merged.{dataset}.stats"), aligner=aln_methods, dataset=["reference", "filtered"])
     output: touch("merged.done")

rule create_list:
     input:
         asms = glob.glob(os.path.join("Assemblies", "{aligner}", "*", "results.gtf"))
     output: os.path.join("CuffMerge", "{aligner}", "gtf_list.txt")
     message: """set +u && mkdir -p $(dirname {output}) && /bin/ls Assemblies/{wildcards.aligner}/*/results.gtf > {output}"""
     shell: """set +u && mkdir -p $(dirname {output}) && cd $(dirname {output}) && /bin/ls ../../Assemblies/{wildcards.aligner}/*/results.gtf > $(basename {output})"""

rule merge:
    input:
         gtf_list=os.path.join("CuffMerge", "{aligner}", "gtf_list.txt"),
	 genome=os.path.join("Reference", "genome.fa")
    output: os.path.join("CuffMerge", "{aligner}", "merged.gtf")
    shell: """set +u && source cufflinks-2.1.1 && cd $(dirname {output}) && cuffmerge -o . -p 10 -s ../../{input.genome} $(basename {input.gtf_list}) && set -u"""

rule compare:
    input:
         merged=os.path.join("CuffMerge", "{aligner}", "merged.gtf"),
	 ref_gff=os.path.join("Reference", "{dataset}.gff3")
    output: os.path.join("CuffMerge", "{aligner}", "Merged.{dataset}.stats")
    log: os.path.join("CuffMerge", "{aligner}", "compare.{dataset}.log")
    shell: """set +u && source mikado-1.0.0b3 && mikado compare -r {input.ref_gff} -p {input.merged} -l {log} -o $(dirname {log})/Merged.{wildcards.dataset} && set -u"""
