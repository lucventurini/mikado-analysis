import os
import glob

aligners = ["STAR", "TopHat"]
methods = ["nosplit", "split", "lenient", "stringent", "permissive"]

rule all:
     input:
         single=expand(os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.stats"), method=methods, aligner=aligners),
         single_filter=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.stats"), method=methods, aligner=aligners),
	 single_standard=expand(os.path.join("Comparisons", "{aligner}", "Original", "Mikado_standard", "{method}.stats"), method=methods, aligner=aligners),
	 single_filter_standard=expand(os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado_standard", "{method}.stats"), method=methods, aligner=aligners)
     output: touch("mikado_compare.done")

rule orig_compare:
     input:
        gtf=os.path.join("Mikado", "{aligner}", "Daijin", "5-mikado", "pick", "{method}", "mikado-{method}.loci.gff3"),
        midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Mikado", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Original", "Mikado"),
        reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && source mikado-devel && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filter_compare:
     input:
         gtf=os.path.join("Mikado", "{aligner}", "Daijin", "5-mikado", "pick", "{method}", "mikado-{method}.loci.gff3"),
         midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado"),
        reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && source mikado-devel && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule orig_standard_compare:
     input:
        gtf=os.path.join("Mikado", "{aligner}", "Daijin", "5-mikado", "pick_standard", "{method}", "mikado-{method}.loci.gff3"),
        midx=os.path.join("Reference", "reference.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Original", "Mikado_standard", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Original", "Mikado_standard", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Original", "Mikado_standard", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Original", "Mikado_standard"),
        reference=os.path.join("Reference", "reference.gff3")
     shell: """set +u && source mikado-devel && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule filter_standard_compare:
     input:
         gtf=os.path.join("Mikado", "{aligner}", "Daijin", "5-mikado", "pick_standard", "{method}", "mikado-{method}.loci.gff3"),
         midx=os.path.join("Reference", "filtered.gff3.midx")
     output: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado_standard", "{method}.stats")
     log: os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado_standard", "{method}.log")
     params:
        prefix=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado_standard", "{method}"),
	out_folder=os.path.join("Comparisons", "{aligner}", "Filtered", "Mikado_standard"),
        reference=os.path.join("Reference", "filtered.gff3")
     shell: """set +u && source mikado-devel && mkdir -p {params.out_folder} && mikado compare -r {params.reference} -p {input.gtf} -l {log} -o {params.prefix} && set -u"""

rule ref_midx:
     input: os.path.join("Reference", "reference.gff3")
     output:
        midx=os.path.join("Reference", "reference.gff3.midx")
     log: os.path.join("Reference", "index.log")
     shell: """set +u && source mikado-devel && mikado compare -r {input} --index -l {log} && set -u"""

rule filtered_midx:
     input: os.path.join("Reference", "filtered.gff3")
     output:
        midx=os.path.join("Reference", "filtered.gff3.midx")
     log: os.path.join("Reference", "index.filtered.log")
     shell: """set +u && source mikado-devel && mikado compare -r {input} --index -l {log} && set -u"""