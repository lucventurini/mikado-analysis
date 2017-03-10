import os
import glob

aligners = ["STAR", "TopHat"]
methods = ["nosplit", "split", "lenient", "stringent", "permissive"]
externals = ["evigene", "stringtie_merge", "cuffmerge"]

rule all:
     input:
         single=expand(os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.stats"), method=methods, aligner=aligners),
         single_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.stats"), method=methods, aligner=aligners),
	 external=expand(os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{external}.stats"), external=externals, aligner=aligners),
	 external_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{external}.stats"), external=externals, aligner=aligners),
	 prepare=expand(os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "prepared.stats"), aligner=aligners),
	 prepare_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "prepared.stats"), aligner=aligners)
     output: touch("mikado_compare.done")

def external_mapper(wildcards):
    base = os.path.join(os.getcwd(), os.path.basename(os.getcwd()), wildcards.aligner)
    if wildcards.external == "evigene":
        return os.path.join(base, "EviGene", "evigene.gtf")
    elif wildcards.external == "stringtie_merge":
        return os.path.join(base, "StringtieMerge", "stringtie.merged.gtf")
    elif wildcards.external == "cuffmerge":
        return os.path.join(base, "Cuffmerge", "merged.gtf")
    else:
        raise ValueError("Invalid wildcard: {}".format(wildcards.external))

rule external_compare:
     input:
       gtf=external_mapper,
       midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{external}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{external}.log")
     params:
       prefix=os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{external}"),
       out_folder=os.path.join("Comparisons", "{aligner}", "Original", "Mikado"),
       reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""
     
rule external_filtered_compare:  
     input:
       gtf=external_mapper,
       midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{external}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{external}.log")
     params:
       prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{external}"),
       out_folder=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado"),
       reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule orig_prepared_compare:
     input:
       gtf=os.path.join(os.path.basename(os.getcwd()), "{aligner}", "5-mikado", "mikado_prepared.gtf"),
       midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "prepared.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "prepared.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "prepared"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Original", "Mikado"),
        reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""     

rule filtered_prepared_compare:
     input:
       gtf=os.path.join(os.path.basename(os.getcwd()), "{aligner}", "5-mikado", "mikado_prepared.gtf"),
       midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "prepared.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "prepared.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "prepared"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado"),
        reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""     

rule orig_compare:
     input:
        gtf=os.path.join(os.path.basename(os.getcwd()), "{aligner}", "5-mikado", "pick", "{method}", "mikado-{method}.loci.gff3"),
        midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Original", "Mikado"),
        reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filter_compare:
     input:
         gtf=os.path.join(os.path.basename(os.getcwd()), "{aligner}", "5-mikado", "pick", "{method}", "mikado-{method}.loci.gff3"),
         midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado"),
        reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && ml mikado/1.0 && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""
     
rule ref_midx:
     input: os.path.join("Reference", "reference.gff3")
     output:
        midx=os.path.join("Reference", "reference.gff3.midx")
     log: os.path.join("Reference", "index.log")
     shell: """set +u && ml mikado/1.0 && mikado compare -r {input} --index -l {log} && set -u"""

rule filtered_midx:
     input: os.path.join("Reference", "filtered.gff3")
     output:
        midx=os.path.join("Reference", "filtered.gff3.midx")
     log: os.path.join("Reference", "index.filtered.log")
     shell: """set +u && ml mikado/1.0 && mikado compare -r {input} --index -l {log} && set -u"""