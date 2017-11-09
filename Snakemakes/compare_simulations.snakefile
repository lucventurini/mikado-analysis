import os
import glob
import itertools

aligners = ["STAR", "TopHat"]
methods = ["CLASS", "Cufflinks", "Stringtie", "Trinity"]

gtfs = dict()
for aligner in aligners:
    gtfs[aligner] = dict()
    for method in methods:
       gtfs[aligner][method] = os.path.join(os.path.basename(os.getcwd()),
                                            "3-assemblies",
                                            "output",
                                            "{asm}-0-{ali}-combined-0.{suff}".format(asm=method.lower(), ali=aligner.lower(), suff="gtf"))

# if os.path.basename(os.getcwd()) == "Athaliana":
#     ss = True
#     ss_flag = "-s"
# else:
#     ss = False
#     ss_flag = ""


ss = False
ss_flag = ""
if os.path.basename(os.getcwd()) in ("Athaliana", "Hsapiens"):
    var = "var_3"
else:
    var = "var_10"

rule all:
     input:
         single=expand(os.path.join("Comparisons", "{aligner}", "Original", "{method}.stats"), method=methods, aligner=aligners),
         single_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.stats"), method=methods, aligner=aligners),
         full=expand(os.path.join("Comparisons", "{aligner}", "Original", "Combined.stats"), aligner=aligners),
         full_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.stats"), aligner=aligners),
         combined=os.path.join("Comparisons", "Combined", "Original", "Combined.stats"),
	 stats=expand(os.path.join("Reference", "{dataset}.stats"), dataset=["reference", "filtered"]),
         combined_filter=os.path.join("Comparisons", "Combined", "Filtered", "Combined.stats"),
	 self_compare=expand(os.path.join("Comparisons", "Self", "{dataset}.stats"), dataset=["reference", "filtered"]),
	 internal_compare=expand(os.path.join("Comparisons", "Internal", "{dataset}.tmap"), dataset=["reference", "filtered"])
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
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Combined")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filtered_compare_combined:
     input:
        midx=os.path.join("Reference", "filtered.gff3.midx"),
        gtf=os.path.join("Comparisons", "{aligner}", "Combined.gtf")
     output:
        os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Combined.log")
     params:
        reference=os.path.join("Reference", "filtered.gff3"),
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Combined")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule combine_all:
     input: list(itertools.chain([list(gtfs[aligner].values()) for aligner in aligners]))
     output: os.path.join("Comparisons", "Combined", "Combined.gtf")
     params:
        comb_folder=os.path.join("Comparisons", "Combined"),
        temp=os.path.join("Comparisons", "Combined", "gtf.tmp")
     message: "Collapsing input transcripts from both aligners with gffread."
     log: os.path.join("logs", "Combination.log")
     shell: """set +u && source cufflinks-2.1.1 && mkdir -p {params.comb_folder} && cat <(sed 's:gene_id ":gene_id "STAR_:; s:transcript_id ":transcript_id "STAR_:;' Comparisons/STAR/Combined.gtf) <(sed 's:gene_id ":gene_id "TopHat_:; s:transcript_id ":transcript_id "TopHat_:;' Comparisons/TopHat/Combined.gtf)  > {params.temp} && gffread -T -M --cluster-only {params.temp} -o {output} && sed -i 's:gene_id ".*; locus ":gene_id ":' {output} && rm -f {params.temp} && sleep 60 && set -u"""

rule orig_compare_all_combined:
     input:
        midx=os.path.join("Reference", "reference.gff3.midx"),     
        gtf=rules.combine_all.output
     output: 
        os.path.join("Comparisons", "Combined", "Original", "Combined.stats")
     log: os.path.join("Comparisons", "Combined", "Original", "Combined.log")
     params:
        reference=os.path.join("Reference", "reference.gff3"),
        prefix=os.path.join("Comparisons", "Combined", "Original", "Combined")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filtered_compare_all_combined:
     input:
        midx=os.path.join("Reference", "filtered.gff3.midx"),     
        gtf=rules.combine_all.output
     output: 
        os.path.join("Comparisons", "Combined", "Filtered", "Combined.stats")
     log: os.path.join("Comparisons", "Combined", "Filtered", "Combined.log")
     params:
        reference=os.path.join("Reference", "filtered.gff3"),
        prefix=os.path.join("Comparisons", "Combined", "Filtered", "Combined")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

def my_aligner_samples(wildcards):
    return list(gtfs[wildcards.aligner].values())

rule prepare_trinity_star:
     input: "{}.gff".format(os.path.splitext(gtfs["STAR"]["Trinity"])[0])
     output:
       gtf=gtfs["STAR"]["Trinity"],
       tmp1=temp(gtfs["STAR"]["Trinity"] + ".tmp1"),
       tmp2=temp(gtfs["STAR"]["Trinity"] + ".tmp2"),
     log: os.path.join("logs", "prepare_trinity_star.log")
     shell: """set +u && source cufflinks-2.2.1 && gffread --force-exons -T -O -o {output.tmp1} {input} && sed -i 's:transcript_id "\([^"]*\)":transcript_id "\\1"; gene_id "\\1.gene":' {output.tmp1} && add_transcript_feature_to_gtf.py {output.tmp1} {output.tmp2} && gffread -O -T -M --cluster-only -o {output.tmp1} {output.tmp2} && sed 's:gene_id ".*; locus ":gene_id ":' {output.tmp1} > {output.gtf} && sleep 30 && set -u"""

rule prepare_trinity_tophat:
     input: "{}.gff".format(os.path.splitext(gtfs["TopHat"]["Trinity"])[0])
     output:
       gtf=gtfs["TopHat"]["Trinity"],
       tmp1=temp(gtfs["TopHat"]["Trinity"] + ".tmp1"),
       tmp2=temp(gtfs["TopHat"]["Trinity"] + ".tmp2"),
     log: os.path.join("logs", "prepare_trinity_tophat.log")
     shell: """set +u && source cufflinks-2.2.1 && gffread --force-exons -T -O -o {output.tmp1} {input} && sed -i 's:transcript_id "\([^"]*\)":transcript_id "\\1"; gene_id "\\1.gene":' {output.tmp1} && add_transcript_feature_to_gtf.py {output.tmp1} {output.tmp2} && gffread -O -T -M --cluster-only -o {output.tmp1} {output.tmp2} && sed 's:gene_id ".*; locus ":gene_id ":' {output.tmp1} > {output.gtf} && sleep 30 && set -u"""

rule combine:
     input: my_aligner_samples
     output: os.path.join("Comparisons", "{aligner}", "Combined.gtf")
     log: os.path.join("logs", "combine.{aligner}.log")
     params:
         temp=os.path.join("Comparisons", "{aligner}", "gtf.temp"),
         my_gtfs=lambda wildcards: " ".join(my_aligner_samples(wildcards))
     message: "Collapsing all the input transcripts with gffread."       
     shell: """set +u && source cufflinks-2.1.1 && cat {params.my_gtfs} | gffread -O -T -M --cluster-only -o {output} - && sed -i 's:gene_id ".*; locus ":gene_id ":' {output} && rm -f {params.temp} && sleep 60"""

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
        reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filter_compare:
     input:
         gtf=my_sample,
         midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "{method}"),
        reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""   

rule ref_midx:
     input: os.path.join("Reference", "reference.gff3")
     output:
        midx=os.path.join("Reference", "reference.gff3.midx")
     log: os.path.join("Reference", "index.log")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {input} --index -l {log} && set -u"""

rule filtered_midx:
     input: os.path.join("Reference", "filtered.gff3")
     output:
        midx=os.path.join("Reference", "filtered.gff3.midx")
     log: os.path.join("Reference", "index.filtered.log")
     shell: """set +u && ml mikado/1.0.1 && mikado compare -r {input} --index -l {log} && set -u"""

rule create_filter:
     input:
        gff3=os.path.join("Reference", "reference.gff3"),
	list=os.path.join("Filter", "final.reconstructable")
     log: "logs/create_filter.log"
     output: os.path.join("Reference", "filtered.gff3")
     shell: """set +u && ml mikado/1.0.1 && mikado util grep {input.list} {input.gff3} {output} && set -u"""

rule create_list:
     input:
     	asm_recon=os.path.join("Filter", "assemblers.reconstructable"),
	ali_recon=os.path.join("Filter", "aligners.reconstructable")
     output:
         list=os.path.join("Filter", "final.reconstructable")
     log: "logs/create_list.log"
     shell: """set +u && cat {input.ali_recon} {input.asm_recon} | sort -u > {output.list}"""

rule create_asm_list:
     input:
         stat=os.path.join("Comparisons", "Combined", "Original", "Combined.stats")
     output: os.path.join("Filter", "assemblers.reconstructable")
     params:
         refmap=os.path.join("Comparisons", "Combined", "Original", "Combined.refmap")
     log: "logs/create_asm_list.log"
     shell: """set +u && mkdir -p Filter/ && cd Filter/ && if [ ! -s $(basename {params.refmap}) ]; then ln -s ../{params.refmap} .; fi && cat *refmap | awk '$2~"(=|_)"' | cut -f 1,8 | sort -u > assemblers.reconstructable"""

rule compare_self:
     input:
         midx=os.path.join("Reference", "{dataset}.gff3.midx")
     output:
         stat=os.path.join("Comparisons", "Self", "{dataset}.stats")
     params:
         gff3=os.path.join("Reference", "{dataset}.gff3")
     log: "Comparisons/Self/{dataset}.log"
     shell: """set +u && mkdir -p Comparisons/Self/ && ml mikado/1.0.1 && mikado compare -r {params.gff3} --self -o Comparisons/Self/{wildcards.dataset} --log Comparisons/Self/{wildcards.dataset}.log && set -u"""

rule compare_internal:
     input:
         midx=os.path.join("Reference", "{dataset}.gff3.midx")
     output:
         stat=os.path.join("Comparisons", "Internal", "{dataset}.tmap")
     params:
         gff3=os.path.join("Reference", "{dataset}.gff3")
     log: "Comparisons/Internal/{dataset}.log"
     shell: """set +u && mkdir -p Comparisons/Self/ && ml mikado/1.0.1 && mikado compare -r {params.gff3} --internal -o Comparisons/Internal/{wildcards.dataset} --log Comparisons/Internal/{wildcards.dataset}.log && set -u"""

rule create_ali_list:
     input:
         bam=glob.glob(os.path.join("datasets", "*{}*bam".format(var)))[0],
	 introns=os.path.join("datasets", "reads", "sim", var, "portcullis", "2-junc", "portcullis.junctions.gff3"),
	 trimmed=os.path.join("Reference", "trimmed.gff3"),
	 fai=os.path.join("Reference", "genome.fa.fai")
     output: os.path.join("Filter", "aligners.reconstructable")
     log: os.path.join("Filter", "aligners.log")
     params:
         prefix=os.path.join("Filter", "aligners"),
     shell: """set +u && ml gcc/5.2.0 bedtools/2.27.0_beta genometools/1.5.9 python/3.5 numpy scipy sklearn blast biopython mikado/1.0 && mkdir -p Filter && python3 ../../get_filtered_reference.py -b {input.bam} -p {input.introns} -o {params.prefix} -r {input.trimmed} --log {log} -gf {input.fai}"""

# rule create_ali_spec_list:
#     input:
#         bam=os.path.join("Alignment", "{aligner}", "alignment.bam"),
# 	introns=os.path.join("Alignment", "{aligner}", "Portcullis-0.17.0", "junc", "portcullis.introns.gff3"),
# 	trimmed=os.path.join("Reference", "trimmed.gff3"),
# 	fai=(os.path.join("Reference", "genome.fa.fai") if os.path.basename(os.getcwd()) != "Hsapiens" else os.path.join("Reference", "genome.fa.correct_order.fai"))
#     output: os.path.join("Filter", "{aligner}", "reconstructable.reconstructable")
#     log: os.path.join("Filter", "{aligner}", "log.txt")
#     params:
#         prefix=os.path.join("Filter", "{aligner}", "reconstructable"),
# 	ss=ss_flag

rule stats:
    input: os.path.join("Reference", "{dataset}.gff3")
    output: os.path.join("Reference", "{dataset}.stats")
    log: os.path.join("logs/{dataset}.stat.log")
    shell: """set +u && ml mikado/1.0.1 && mikado util stats {input} {output}"""

# rule assembly_stats:
#     input: os.path.join("Assemblies", "{aligner}", "{method}", "results.gtf")
#     output: os.path.join("Assemblies", "{aligner}", "{method}", "results.stats")
#     shell: "set +u && source mikado-1.0.0b3 && mikado util stats {input} {output} && set -u"