import os
import glob
import itertools

aligners = ["STAR", "TopHat"]
methods = ["CLASS", "Cufflinks", "Stringtie", "Trinity"]

methods = ["CLASS", "Cufflinks", "Stringtie", "Trinity"]

def my_sample(wildcards):
    return gtfs[wildcards.aligner][wildcards.method]

def my_aligner_samples(wildcards):
    return list(gtfs[wildcards.aligner].values())

gtfs = dict()
for aligner in aligners:
    gtfs[aligner] = dict()
    for method in methods:
       gtfs[aligner][method] = os.path.join(os.path.basename(os.getcwd()),
                                            "3-assemblies",
                                            "output",
                                            "{asm}-0-{ali}-combined-0.{suff}".format(asm=method.lower(), ali=aligner.lower(), suff="gtf"))

if os.path.basename(os.getcwd()) == "Athaliana":
    ss = True
    ss_flag = "-s"
else:
    ss = False
    ss_flag = ""

rule all:
     input:
         single=expand(os.path.join("Comparisons", "{aligner}", "Original", "{method}.stats"), method=methods, aligner=aligners),
         single_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.stats"), method=methods, aligner=aligners),
         full=expand(os.path.join("Comparisons", "{aligner}", "Original", "Combined.stats"), aligner=aligners),
         full_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.stats"), aligner=aligners),
         combined=os.path.join("Comparisons", "Combined", "Original", "Combined.stats"),
	 stats=expand(os.path.join("Reference", "{dataset}.stats"), dataset=["reference", "filtered"]),
	 assembly_stats=expand(os.path.join("Assemblies", "{aligner}", "{method}", "results.stats"), aligner=aligners, method=methods),
         combined_filter=os.path.join("Comparisons", "Combined", "Filtered", "Combined.stats"),
	 self_compare=expand(os.path.join("Comparisons", "Self", "{dataset}.stats"), dataset=["reference", "filtered"]),
	 internal_compare=expand(os.path.join("Comparisons", "Internal", "{dataset}.tmap"), dataset=["reference", "filtered"])
     params:
       memory=2000
     output: touch("all.done")


rule orig_compare_combined:
     input:
        midx=os.path.join("Reference", "reference.gff3.midx"),
        gtf=os.path.join("Comparisons", "{aligner}", "Combined.gtf")
     output:
        os.path.join("Comparisons", "{aligner}", "Original", "Combined.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Combined.log")
     params:
        reference=os.path.join("Reference", "reference.gff3"),
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Combined"),
	memory=15000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filtered_compare_combined:
     input:
        midx=os.path.join("Reference", "filtered.gff3.midx"),
        gtf=os.path.join("Comparisons", "{aligner}", "Combined.gtf")
     output:
        os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.log")
     params:
        reference=os.path.join("Reference", "filtered.gff3"),
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Combined"),
	memory=15000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule combine_all:
     input:
       star_gtf=os.path.join("Comparisons", "STAR", "Combined.gtf"),
       top_gtf=os.path.join("Comparisons", "TopHat", "Combined.gtf")
     output:
       gtf=os.path.join("Comparisons", "Combined", "Combined.gtf"),
       tmp=temp(os.path.join("Comparisons", "Combined", "tmp.gtf"))
     params:
        comb_folder=os.path.join("Comparisons", "Combined"),
	memory=10000
     shell: """set +u && ml Cufflinks/2.2.2-foss-2016a && mkdir -p {params.comb_folder} && cat <(sed 's:gene_id ":gene_id "STAR_:; s:transcript_id ":transcript_id "STAR_:;' {input.star_gtf}) <(sed 's:gene_id ":gene_id "TopHat_:; s:transcript_id ":transcript_id "TopHat_:;' {input.top_gtf})  > {output.tmp} && gffread -T -M --cluster-only -o {output.gtf} {output.tmp} && sed -i 's:gene_id ".*; locus ":gene_id ":' {output.gtf} && sleep 60 && set -u"""

rule orig_compare_all_combined:
     input:
        midx=os.path.join("Reference", "reference.gff3.midx"),     
        gtf=rules.combine_all.output.gtf
     output: 
        os.path.join("Comparisons", "Combined", "Original", "Combined.stats")
     log: os.path.join("Comparisons", "Combined", "Original", "Combined.log")
     params:
        reference=os.path.join("Reference", "reference.gff3"),
        prefix=os.path.join("Comparisons", "Combined", "Original", "Combined"),
        memory=20000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filtered_compare_all_combined:
     input:
        midx=os.path.join("Reference", "filtered.gff3.midx"),     
        gtf=rules.combine_all.output.gtf
     output: 
        os.path.join("Comparisons", "Combined", "Filtered", "Combined.stats")
     log: os.path.join("Comparisons", "Combined", "Filtered", "Combined.log")
     params:
        reference=os.path.join("Reference", "filtered.gff3"),
        prefix=os.path.join("Comparisons", "Combined", "Filtered", "Combined"),
        memory=20000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

def my_aligner_samples(wildcards):
    return list(gtfs[wildcards.aligner].values())

rule prepare_trinity_star:
     input: "{}.gff".format(os.path.splitext(gtfs["STAR"]["Trinity"])[0])
     output:
       gtf=gtfs["STAR"]["Trinity"],
       tmp1=temp(gtfs["STAR"]["Trinity"] + ".tmp1"),
       tmp2=temp(gtfs["STAR"]["Trinity"] + ".tmp2"),
     params:
       memory=7000
     shell: """set +u && ml Cufflinks/2.2.2-foss-2016a && gffread --force-exons -T -O -o {output.tmp1} {input} && sed -i 's:transcript_id "\([^"]*\)":transcript_id "\\1"; gene_id "\\1.gene":' {output.tmp1} && add_transcript_feature_to_gtf.py {output.tmp1} {output.tmp2} && gffread -O -T -M --cluster-only -o {output.tmp1} {output.tmp2} && sed 's:gene_id ".*; locus ":gene_id ":' {output.tmp1} > {output.gtf} && sleep 30 && set -u"""

rule prepare_trinity_tophat:
     input: "{}.gff".format(os.path.splitext(gtfs["TopHat"]["Trinity"])[0])
     output:
       gtf=gtfs["TopHat"]["Trinity"],
       tmp1=temp(gtfs["TopHat"]["Trinity"] + ".tmp1"),
       tmp2=temp(gtfs["TopHat"]["Trinity"] + ".tmp2"),
     params:
       memory=7000
     shell: """set +u && ml Cufflinks/2.2.2-foss-2016a && gffread --force-exons -T -O -o {output.tmp1} {input} && sed -i 's:transcript_id "\([^"]*\)":transcript_id "\\1"; gene_id "\\1.gene":' {output.tmp1} && add_transcript_feature_to_gtf.py {output.tmp1} {output.tmp2} && gffread -O -T -M --cluster-only -o {output.tmp1} {output.tmp2} && sed 's:gene_id ".*; locus ":gene_id ":' {output.tmp1} > {output.gtf} && sleep 30 && set -u"""

rule combine:
     input: my_aligner_samples
     output: os.path.join("Comparisons", "{aligner}", "Combined.gtf")
     params:
         temp=os.path.join("Comparisons", "{aligner}", "gtf.temp"),
         my_gtfs=lambda wildcards: " ".join(my_aligner_samples(wildcards)),
         memory=10000
     message: "Collapsing all the input transcripts with gffread."
     shell: """set +u && ml Cufflinks/2.2.2-foss-2016a && cat {params.my_gtfs} | gffread -O -T -M --cluster-only -o {output} - && sed -i 's:gene_id ".*; locus ":gene_id ":' {output} && rm -f {params.temp} && sleep 60"""

def my_sample(wildcards):
    return gtfs[wildcards.aligner][wildcards.method]

rule orig_compare:
     input:
         gtf=my_sample,
         midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "{method}"),
        reference=os.path.join("Reference", "reference.gff3"),
	memory=15000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filter_compare:
     input:
         gtf=my_sample,
         midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "{method}"),
        reference=os.path.join("Reference", "filtered.gff3"),
	memory=15000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""   

rule ref_midx:
     input: os.path.join("Reference", "reference.gff3")
     output:
        midx=os.path.join("Reference", "reference.gff3.midx")
     log: os.path.join("Reference", "index.log")
     params:
       memory=20000
     shell: """set +u && ml mikado/1.0 && mikado compare -r {input} --index -l {log} && set -u"""

rule filtered_midx:
     input: os.path.join("Reference", "filtered.gff3")
     output:
        midx=os.path.join("Reference", "filtered.gff3.midx")
     params:
        memory=15000
     log: os.path.join("Reference", "index.filtered.log")
     shell: """set +u && ml mikado/1.0 && mikado compare -r {input} --index -l {log} && set -u"""

rule create_filter:
     input:
        gff3=os.path.join("Reference", "reference.gff3"),
	list=os.path.join("Filter", "final.reconstructable")
     output: os.path.join("Reference", "filtered.gff3")
     params:
        memory=3000
     shell: """set +u && ml mikado/1.0 && mikado util grep {input.list} {input.gff3} {output} && set -u"""

rule create_list:
     input:
     	asm_recon=os.path.join("Filter", "assemblers.reconstructable"),
	ali_recon=os.path.join("Filter", "aligners.reconstructable")
     output:
         list=os.path.join("Filter", "final.reconstructable")
     params:
         memory=5000
     shell: """set +u && cat {input.ali_recon} {input.asm_recon} | sort -u > {output.list}"""

rule create_asm_list:
     input:
         stat=os.path.join("Comparisons", "Combined", "Original", "Combined.stats")
     output: os.path.join("Filter", "assemblers.reconstructable")
     params:
         refmap=os.path.join("Comparisons", "Combined", "Original", "Combined.refmap"),
	 memory=3000
     shell: """set +u && mkdir -p Filter/ && cd Filter/ && if [ ! -s $(basename {params.refmap}) ]; then ln -s ../{params.refmap} .; fi && cat *refmap | awk '$2~"(=|_)"' | cut -f 1,8 | sort -u > assemblers.reconstructable"""

rule create_ali_list:
     input:
         star_list=os.path.join("Filter", "STAR", "reconstructable.reconstructable"),
	 top_list=os.path.join("Filter", "TopHat", "reconstructable.reconstructable")
     output: os.path.join("Filter", "aligners.reconstructable")
     params:
       memory=5000
     shell: """set +u && cat {input.star_list} {input.top_list} | sort -u > {output}"""

def to_bam(wildcards):
    return os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "2-alignments", "output", "{}-combined-0.sorted.bam".format(wildcards.aligner.lower()))

def to_portcullis_gff3(wildcards):
    return os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "4-portcullis", "portcullis_{}-combined-0".format(wildcards.aligner.lower()), "2-junc", "portcullis.junctions.gff3")

def to_portcullis_tab(wildcards):
    return os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "4-portcullis", "portcullis_{}-combined-0".format(wildcards.aligner.lower()), "2-junc", "portcullis.junctions.tab")

rule convert_portcullis_tab_to_gff3:
    input: to_portcullis_tab
    output: os.path.join(os.getcwd(), os.path.basename(os.getcwd()), "4-portcullis", "portcullis_{}-combined-0".format("{aligner}".lower()), "2-junc", "portcullis.junctions.gff3")
    params:
      memory=7000
    shell: """set +u && ml portcullis/1.0.0_beta5 && junctools convert -if portcullis -of igff -o {output} {input}"""

rule create_ali_spec_list:
    input:
        bam=to_bam,
	introns=to_portcullis_gff3,
	trimmed=os.path.join("Reference", "trimmed.gff3"),
	fai=(os.path.join("Reference", "genome.fa.fai") if os.path.basename(os.getcwd()) != "Hsapiens" else os.path.join("Reference", "genome.fa.correct_order.fai"))
    output: os.path.join("Filter", "{aligner}", "reconstructable.reconstructable")
    log: os.path.join("Filter", "{aligner}", "log.txt")
    params:
        prefix=os.path.join("Filter", "{aligner}", "reconstructable"),
	ss=ss_flag,
	memory=100000
    shell: """set +u && ml gcc/5.2.0 bedtools/2.27.0_beta genometools/1.5.9 python/3.5 numpy scipy sklearn blast biopython mikado/1.0 && mkdir -p $(dirname {output}) && python3 ../get_filtered_reference.py -b {input.bam} -p {input.introns} -o {params.prefix} -r {input.trimmed} {params.ss} --log {log} -gf {input.fai}"""

rule compare_self:
     input:
         midx=os.path.join("Reference", "{dataset}.gff3.midx")
     output:
         stat=os.path.join("Comparisons", "Self", "{dataset}.stats")
     params:
         gff3=os.path.join("Reference", "{dataset}.gff3"),
	 memory=20000
     shell: """set +u && mkdir -p Comparisons/Self/ && ml mikado/1.0 && mikado compare -r {params.gff3} --self -o Comparisons/Self/{wildcards.dataset} --log Comparisons/Self/{wildcards.dataset}.log && set -u"""

rule compare_internal:
     input:
         midx=os.path.join("Reference", "{dataset}.gff3.midx")
     output:
         stat=os.path.join("Comparisons", "Internal", "{dataset}.tmap")
     params:
         gff3=os.path.join("Reference", "{dataset}.gff3"),
	 memory=10000
     shell: """set +u && mkdir -p Comparisons/Self/ && ml mikado/1.0 && mikado compare -r {params.gff3} --internal -o Comparisons/Internal/{wildcards.dataset} --log Comparisons/Internal/{wildcards.dataset}.log && set -u"""

rule stats:
    input: os.path.join("Reference", "{dataset}.gff3")
    output: os.path.join("Reference", "{dataset}.stats")
    params:
      memory=10000
    shell: """set +u && ml mikado/1.0 && mikado util stats {input} {output}"""
