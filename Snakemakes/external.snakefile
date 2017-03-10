import os
import sys


aligners = ["STAR", "TopHat"]
# methods = ["nosplit", "split", "lenient", "stringent", "permissive"]

rule all:
     input:
     	 stringtie_merge=expand(os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "StringtieMerge", "stringtie.merged.gtf"), aligner=aligners),
         evigene=expand(os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "evigene.gtf"), aligner=aligners),
	 cuffmerge=expand(os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "Cuffmerge", "merged.gtf"), aligner=aligners),
     output: touch("external.done")

rule cuffmerge:
  input:
     gtf=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.gtf"),
     genome=os.path.join(os.getcwd(), "Reference", "genome.fa")
  output: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "Cuffmerge", "merged.gtf")
  log: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}",  "Cuffmerge", "cuffmerge.log")
  threads: 10
  params:
    out_folder=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "Cuffmerge")
  shell: """set +u && mkdir -p {params.out_folder} && cd {params.out_folder} && ml Cufflinks/2.2.2-foss-2016a python/2.7 && cuffmerge -o . -s {input.genome} -p {threads} <(echo {input.gtf}) 2>&1 > cuffmerge.log"""

rule gtf_to_fasta:
  input:
    gtf=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.gtf"),
    genome=os.path.join(os.getcwd(), "Reference", "genome.fa")
  output:
    fasta=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.fasta")
  shell: """set +u && ml Cufflinks/2.2.2-foss-2016a SAMtools/1.3-foss-2016a && gffread -g {input.genome} -w {output.fasta} {input.gtf} && samtools faidx {output.fasta} && set -u"""

rule evigene:
  input:
    fasta=rules.gtf_to_fasta.output.fasta
  output:
    okalt=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okalt.fasta"),
    okay=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okay.fasta")
  log: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "evigene.log")
  params:
    out_folder=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene")
  threads: 10
  shell: """set +u && mkdir -p {params.out_folder} && cd {params.out_folder} && source evidentialGene-20160320 && ml gcc/4.9.1 CD-HIT/4.6.4-foss-2016a-2015-0603 && tr2aacds -NCPU={threads} -MAXMEM=7000 -logfile {log} -tidyup -mrnaseq {input.fasta} && set -u"""

rule evigene_list:
  input: gtf=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.gtf")
  output: os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.list.txt")
  shell: """mkdir -p $(dirname ${output}) && cat {input} | sed 's/.*transcript_id "\([^;]*\)"; gene_id "\([^;]*\)".*/\\1\\t\\2/; s/.*gene_id "\([^;]*\)"; transcript_id "\([^;]*\)".*/\\2\\t\\1/;' | uniq | sort -u > {output}"""

rule evigene_gtf:
  input:
    okalt=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okalt.fasta"),
    okay=os.path.join(os.getcwd(),os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okay.fasta"),
    list=rules.evigene_list.output,
    gtf=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.gtf")
  output:
    okalt_fai=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okalt.fasta.fai"),
    okay_fai=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "okayset", "Combined.okay.fasta.fai"),
    tmp=temp(os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "evigene.gtf.tmp")),
    gtf=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "EviGene", "evigene.gtf")
  params:
    out_folder=os.path.join(os.path.basename(os.getcwd()), "{aligner}", "EviGene")
  shell: """set +u && ml samtools && ml python/3.5 biopython blast numpy scipy sklearn mikado && samtools faidx {input.okay} && samtools faidx {input.okalt} && cd {params.out_folder} && mikado util grep <(grep.py <(cut -f 1 {output.okalt_fai} {output.okay_fai}) {input.list}) {input.gtf} {output.tmp} && ml Cufflinks/2.2.2-foss-2016a python/2.7 && gffread -T -o {output.gtf} -M --cluster-only {output.tmp} && sed -i 's:gene_id "[^"]*"; locus :gene_id :' {output.gtf} && set -u"""

rule stringtie_merge:
  input:
    gtf=os.path.join(os.getcwd(), "Comparisons", "{aligner}", "Combined.gtf")
  output: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "StringtieMerge", "stringtie.merged.gtf")
  threads: 1
  log: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "StringtieMerge", "stringtie_merge.log")
  params:
    out_folder=os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "{aligner}", "StringtieMerge")
  shell: """set +u && mkdir -p {params.out_folder} && cd {params.out_folder} && ml stringtie && stringtie -f 0.05 --merge -o {output} {input.gtf} 2>&1 > {log} && set -u"""
