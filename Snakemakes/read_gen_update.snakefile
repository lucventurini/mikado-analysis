import sys

from os import listdir
from os.path import isfile, join, abspath

from snakemake.utils import min_version

min_version("3.4.2")

R1 = config["r1"]
R2 = config["r2"]
R1_NAME, ext1 = os.path.splitext(os.path.basename(R1))
R2_NAME, ext2 = os.path.splitext(os.path.basename(R2))
REF = config["ref_fa"]
REF_GTF = config["ref_gtf"]
NAME = config["name"]
OUT_DIR = config["out_dir"]
OUT_DIR_FULL = os.path.abspath(config["out_dir"])
MIN_INTRON = config["min_intron"]
MAX_INTRON = config["max_intron"]
STRANDEDNESS = config["strandedness"]
THREADS = int(config["threads"])
FRAG_LENGTH = config["frag_length"]


INDEX_STAR_EXTRA = config["index_star_extra"]
ALIGN_TOPHAT_EXTRA = config["align_tophat_extra"]
ALIGN_GSNAP_EXTRA = config["align_gsnap_extra"]
ALIGN_STAR_EXTRA = config["align_star_extra"]
ASM_CUFFLINKS_EXTRA = config["asm_cufflinks_extra"]
ASM_STRINGTIE_EXTRA = config["asm_stringtie_extra"]
JUNC_PORTCULLIS_PREP_EXTRA = config["junc_portcullis_prep_extra"]
JUNC_PORTCULLIS_JUNC_EXTRA = config["junc_portcullis_junc_extra"]


INPUT_SETS = ["real","sim_fixed","sim_var"]
ALIGNMENT_METHODS = ["tophat","star","gsnap","hisat"]
ALIGNMENT_WREF_METHODS = ["tophat","star"]
ASSEMBLY_METHODS = config["asm_methods"]

#Shortcuts
READS_DIR = OUT_DIR + "/reads"
READS_DIR_FULL = os.path.abspath(READS_DIR)
TRUE_DIR = OUT_DIR + "/true"
SIM_DIR = READS_DIR + "/sim"
SIM_DIR_FULL = os.path.abspath(SIM_DIR)
ALIGN_DIR = OUT_DIR + "/alignments"
ALIGN_DIR_FULL = os.path.abspath(ALIGN_DIR)
REF_ALIGN_DIR = READS_DIR + "/alignments"
REF_ALIGN_DIR_FULL = os.path.abspath(REF_ALIGN_DIR)
REF_ASM_DIR = READS_DIR + "/assemblies"
REF_ASM_DIR_FULL = os.path.abspath(REF_ASM_DIR)

CWD = os.getcwd()


CHUNKS = []
for i in range(int(config['chunks'])):
	CHUNKS.append(format(i, '03d'))


HISAT_STRAND = "--rna-strandness=RF" if STRANDEDNESS == "fr-firststrand" else "--rna-strandness=FR" if STRANDEDNESS == "fr-secondstrand" else "F"
PORTCULLIS_STRAND = "firststrand" if STRANDEDNESS == "fr-firststrand" else "secondstrand" if STRANDEDNESS == "fr-secondstrand" else "unstranded"

# Get the last read length in the array, assume the array was sorted and the last element is the largest
MAX_RL = config['read_length'][-1]

RUNS = []
for d in config['depth']:
	srun = "sim_d"+str(d)+"_r"+str(MAX_RL)
	RUNS.append(srun)

for r in config['read_length'][:-1]:
	srun = "sim_d1_r"+str(r)
	RUNS.append(srun)

def getDepthFromRun(srun):
	if srun.startswith("sim"):
		return srun.split("_")[1][1:]
	else:
		return "1"

def getReadLenFromRun(srun):
	if srun.startswith("sim"):
		return srun.split("_")[2][1:]
	else:
		return MAX_RL


def loadPre(command): 
        cc = command.strip() 
        if not cc: 
                return "" 
        else: 
                return "set +u && {} &&".format(cc)


RUN_TIME = config["run_time"]





#########################
Rules

localrules: all, tidy, asm_cufflinks_cov, depth_multipler, hisat_getref_ss, create_chunks

# Define 
rule all:
	input:
		expand([OUT_DIR+"/"+NAME+".{srun}.bam", 
			OUT_DIR+"/"+NAME+".{srun}.bed",
			OUT_DIR+"/"+NAME+".{srun}.tab",
			OUT_DIR+"/"+NAME+".{srun}.R1.fq.gz",
			OUT_DIR+"/"+NAME+".{srun}.R2.fq.gz",
			SIM_DIR+"/{srun}/fastqc/"+NAME+".{srun}.R1_fastqc.html"], srun=RUNS),
		READS_DIR+"/trim/fastqc/"+NAME+".real.R1_fastqc.html",
		REF_ALIGN_DIR+"/hisat/hisat.sorted.bam.bai"


		

rule trim_reads:
	input: 
		r1=R1,
		r2=R2
	output:
		linkr1=OUT_DIR+"/"+NAME+".real.R1.fq.gz",
		linkr2=OUT_DIR+"/"+NAME+".real.R2.fq.gz"
	params:
		load=loadPre(config['load']['trimgalore']),
		outdir=READS_DIR+"/trim",
		r1="reads/trim/"+R1_NAME+"_val_1.fq.gz",
		r2="reads/trim/"+R2_NAME+"_val_2.fq.gz",
		readlen=lambda wildcards: getReadLenFromRun(srun)
	log: READS_DIR+"/trim/trim.log"		
	shell: "{params.load} trim_galore --paired --length {params.readlen} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"


rule fastqc_in:
	input:
		r1=rules.trim_reads.output.linkr1,
		r2=rules.trim_reads.output.linkr2
	output: 
		r1=READS_DIR+"/trim/fastqc/"+NAME+".real.R1_fastqc.html",
		r2=READS_DIR+"/trim/fastqc/"+NAME+".real.R2_fastqc.html"
	params:
		load=loadPre(config['load']['fastqc']),
		outdir=READS_DIR+"/trim/fastqc"
	log: READS_DIR+"/trim/fastqc/fastqc.log"
	shell: "{params.load} fastqc -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1"

rule fastqc_out:
	input:
		r1=SIM_DIR+"/{srun}/spanki/"+NAME+".{srun}.R1.fq.gz",
		r2=SIM_DIR+"/{srun}/spanki/"+NAME+".{srun}.R2.fq.gz"
	output: 
		r1=SIM_DIR+"/{srun}/fastqc/"+NAME+".{srun}.R1_fastqc.html",
		r2=SIM_DIR+"/{srun}/fastqc/"+NAME+".{srun}.R2_fastqc.html"
	params:
		load=loadPre(config['load']['fastqc']),
		outdir=SIM_DIR+"/{srun}/fastqc"
	log: SIM_DIR+"/{srun}/fastqc/fastqc.{srun}.log"
	shell: "{params.load} fastqc -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1"

rule align_bowtie_index:
	input:
		REF
	output: 
		REF_ALIGN_DIR+"/bowtie/index/"+REF+".1.ebwt"
	params: 
		load=loadPre(config['load']['bowtie']),
		outdir=REF_ALIGN_DIR+"/bowtie/index/"+REF
	log:
		REF_ALIGN_DIR+"/bowtie.index.log"
	threads: 1
	message: "Creating bowtie index of genome"
	shell: "{params.load} {RUN_TIME} bowtie-build {input} {params.outdir} > {log} 2>&1"


rule align_bowtie:
	input:
		r1=rules.trim_reads.output.linkr1,
                r2=rules.trim_reads.output.linkr2,
                index=rules.align_bowtie_index.output
	output:
		map=REF_ALIGN_DIR+"/bowtie/bowtie.map"
	params: 
		load=loadPre(config['load']['bowtie']),
		indexdir=REF_ALIGN_DIR+"/bowtie/index/"+REF
	log: REF_ALIGN_DIR+"/bowtie.align.log"		
	threads: THREADS
	shell: "{params.load} {RUN_TIME} bowtie -p {threads} --chunkmbs 1000 -X 500 {params.indexdir} -1 <(zcat {input.r1}) -2 <(zcat {input.r2}) > {output.map} 2> {log}"
	



rule align_hisat_index:
        input: REF
        output: REF_ALIGN_DIR+"/hisat/index/"+NAME+".4.ht2"
        params: load=loadPre(config["load"]["hisat"])
        log: REF_ALIGN_DIR+"/hisat.index.log"
        threads: 1
        message: "Indexing genome with hisat"
        shell: "{params.load} {RUN_TIME} hisat2-build {input} {REF_ALIGN_DIR}/hisat/index/{NAME} > {log} 2>&1"




rule hisat_getref_ss:
	input: REF_GTF
	output: REF_ALIGN_DIR+"/hisat/ref_splicesites.txt"
	params:
                load=loadPre(config["load"]["hisat"]), 
	shell: "{params.load} hisat2_extract_splice_sites.py {input} > {output}"


rule align_hisat: 
        input: 
                r1=rules.trim_reads.output.linkr1, 
                r2=rules.trim_reads.output.linkr2,
		ss=rules.hisat_getref_ss.output,
                index=rules.align_hisat_index.output 
        output: REF_ALIGN_DIR+"/hisat/hisat.unsorted.bam", 
        params: 
                load_h=loadPre(config["load"]["hisat"]), 
                load_s=loadPre(config["load"]["samtools"]), 
                indexdir=REF_ALIGN_DIR+"/hisat/index/"+NAME,
                strand=HISAT_STRAND,
		sam=REF_ALIGN_DIR+"/hisat/hisat.unsorted.sam", 
        log: REF_ALIGN_DIR+"/hisat.align.log" 
        threads: THREADS 
        message: "Aligning input with hisat"
        shell: "{params.load_h} {params.load_s} {RUN_TIME} hisat2 -p {threads} --known-splicesite-infile={input.ss} --dta-cufflinks --min-intronlen={MIN_INTRON} --max-intronlen={MAX_INTRON} {params.strand} -x {params.indexdir} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -@ {threads} -bS > {output}"


rule bam_sort_wref:
	input: rules.align_hisat.output
	output: REF_ALIGN_DIR+"/hisat/hisat.sorted.bam"
	params:
		load=loadPre(config['load']['samtools']),
		temp=REF_ALIGN_DIR+"/hisat/sort_chunk"
	threads: THREADS
	log: REF_ALIGN_DIR+"/hisat/sort.log"
	message: "Using samtools to sort {input}"
	shell: "{params.load} {RUN_TIME} samtools sort -o {output} -O bam -m 2G -T {params.temp} -@ {threads} {input} > {log} 2>&1"

rule bam_index_wref:
	input: rules.bam_sort_wref.output
	output: REF_ALIGN_DIR+"/hisat/hisat.sorted.bam.bai"
	params:
		load=loadPre(config['load']['samtools'])
	threads: 1
	log: REF_ALIGN_DIR+"/hisat/index.log"
	message: "Using samtools to index: {input}"
	shell: "{params.load} {RUN_TIME} samtools index {input} > {log} 2>&1"

rule asm_cufflinks_wref:
	input:
		bam=rules.bam_sort_wref.output,
		ref=REF_GTF
	output:	REF_ASM_DIR+"/cufflinks-hisat/isoforms.fpkm_tracking"
	params:	
		load=loadPre(config['load']['cufflinks']),
		outdir=REF_ASM_DIR+"/cufflinks-hisat"
	log: REF_ASM_DIR+"/cufflinks-hisat.log"
	threads: THREADS
	message: "Using cufflinks to assemble: {input.bam}"
	shell: "{params.load} {RUN_TIME} cufflinks --output-dir={params.outdir} --num-threads={threads} --library-type={STRANDEDNESS} --min-intron-length={MIN_INTRON} --max-intron-length={MAX_INTRON} -G {input.ref} -g {input.ref} --no-update-check {ASM_CUFFLINKS_EXTRA} {input.bam} > {log} 2>&1"

rule asm_cufflinks_cov:
	input: rules.asm_cufflinks_wref.output
	output: REF_ASM_DIR+"/cufflinks-hisat/isoforms.fpkm_tracking_inref"
	threads: 1
	message: "Creating transcript coverage output file"
	shell: "grep -v ^CUFF {input} > {output}"


rule depth_multipler:
	input: rules.asm_cufflinks_cov.output
	output: REF_ASM_DIR+"/cufflinks-hisat/transcripts.{srun}.cov"
	params:	
		load=loadPre(config['load']['portcullis_scripts']),
		depth=lambda wildcards: getDepthFromRun(wildcards.srun)
	message: "Multiplying expression levels as requested"
	shell: "{params.load} cufflinks2spankicov.py -m {params.depth} {input} > {output}"

rule simgen_model:
	input:
		map=rules.align_bowtie.output.map
	output: SIM_DIR+"/model/logfile.txt"
	params: 
		load=loadPre(config['load']['spanki']),
		outdir=SIM_DIR+"/model"
	log: SIM_DIR+"/spanki_model.log"
	threads: 1
	shell: "{params.load} spankisim_models -i {input.map} -e 2 -l {MAX_RL} -o {params.outdir} > {log} 2>&1"


rule sim_fixed_reads:
	input:
		gtf=REF_GTF,
                fa=REF,
		model=rules.simgen_model.output
	output: 
		bam=SIM_DIR_FULL+"/fixed/sim.bam",
		linkr1=READS_DIR+"/sim_fixed.R1.fq",
		linkr2=READS_DIR+"/sim_fixed.R2.fq"
	params: 
		load=loadPre(config['load']['spanki']),
		load_cuff=loadPre(config['load']['cufflinks']),
		outdir=SIM_DIR+"/fixed",
		mdir=SIM_DIR+"/model",
		r1=SIM_DIR_FULL+"/fixed/sim_1.fastq",
		r2=SIM_DIR_FULL+"/fixed/sim_2.fastq",
		readlen=lambda wildcards: getReadLenFromRun(srun)
	log: SIM_DIR+"/spanki_readgen_fixed.log"
	threads: 1	
	shell: "{params.load} {params.load_cuff} spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -cov 30 -bp {params.readlen} -ends 2 -frag 250 > {log} 2>&1 && ln -sf {params.r1} {output.linkr1} && ln -sf {params.r2} {output.linkr2} && touch -h {output.linkr1} {output.linkr2}"
	

rule create_chunks:
        input: cov=rules.depth_multipler.output
        output: SIM_DIR+"/{srun}/spanki/chunks/chunking.done"
        params: 
                load=loadPre(config['load']['portcullis_scripts']),
                nb_chunks=config['chunks'],
                prefix=SIM_DIR+"/{srun}/spanki/chunks/chunk_"
        log: SIM_DIR+"/{srun}/spanki/chunks/chunking.log"
        shell: "{params.load} split_spankicov.py -n {params.nb_chunks} -o {params.prefix} {input} > {log} 2>&1 && touch {output}"



rule sim_var_chunk:
	input:
		gtf=REF_GTF,
                fa=REF,
                transcript_cov=rules.create_chunks.output,
		model=rules.simgen_model.output
	output:
		bam=temp(SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}/sim.bam"),
		r1=temp(SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}/sim_1.fastq.gz"),
		r2=temp(SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}/sim_2.fastq.gz")
	params: 
		load=loadPre(config['load']['spanki']),
		load_cuff=loadPre(config['load']['cufflinks']),
		outdir=SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}",
		cov=SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}.cov",
		sam=SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}/sim.sam",
		mdir=SIM_DIR+"/model",
		readlen=lambda wildcards: getReadLenFromRun(wildcards.srun)
	log: SIM_DIR+"/{srun}/spanki/chunks/chunk_{chunk}/spanki_readgen_var.log"
	threads: 1
	shell: "{params.load} {params.load_cuff} {RUN_TIME} spankisim_transcripts -g {input.gtf} -f {input.fa} -o {params.outdir} -m custom -mdir {params.mdir} -t {params.cov} -bp {params.readlen} -ends 2 -frag {FRAG_LENGTH} > {log} 2>&1 && rm -f {params.sam} && gzip {params.outdir}/sim_?.fastq"


rule merge_fq1:
	input: expand(SIM_DIR+"/{{srun}}/spanki/chunks/chunk_{chunk}/sim_1.fastq.gz", chunk=CHUNKS)
	output: SIM_DIR+"/{srun}/spanki/"+NAME+".{srun}.R1.fq.gz"
	threads: 1
	message: "Merging R1"
	shell: "cat {input} > {output}"

rule merge_fq2:
     input: expand(SIM_DIR+"/{{srun}}/spanki/chunks/chunk_{chunk}/sim_2.fastq.gz", chunk=CHUNKS)
     output: SIM_DIR+"/{srun}/spanki/"+NAME+".{srun}.R2.fq.gz"
     threads: 1
     message: "Merging R2"
     shell: "cat {input} > {output}"

rule merge_bam:
	input: expand(SIM_DIR+"/{{srun}}/spanki/chunks/chunk_{chunk}/sim.bam", chunk=CHUNKS)
	output: SIM_DIR+"/{srun}/spanki/"+NAME+".{srun}.bam"
	params:
		load=loadPre(config['load']['samtools'])
	threads: THREADS
	shell: "{params.load} samtools merge -@ {threads} -f {output} {input}"


	


rule portcullis_prep:
	input:
		bam=rules.merge_bam.output,
		ref=REF
	output: 
		bai=SIM_DIR+"/{srun}/portcullis/1-prep/portcullis.sorted.alignments.bam.bai"
	params:	
		load=loadPre(config['load']['portcullis']),
		outdir=SIM_DIR+"/{srun}/portcullis/1-prep"
	log: SIM_DIR+"/{srun}/portcullis/portcullis.{srun}.1-prep.log"
	threads: THREADS
	message: "Preparing portcullis to extract junctions from: {input}"
	shell: "{params.load} {RUN_TIME} portcullis prep -o {params.outdir} -t {threads} {input.ref} {input.bam} > {log} 2>&1"

rule portcullis_junc:
	input:
		rules.portcullis_prep.output
	output: 
		bed=SIM_DIR+"/{srun}/portcullis/2-junc/portcullis.junctions.bed",
		tab=SIM_DIR+"/{srun}/portcullis/2-junc/portcullis.junctions.tab"
	params:	
		load=loadPre(config['load']['portcullis']),
		outdir=SIM_DIR+"/{srun}/portcullis/2-junc",
		prepdir=SIM_DIR+"/{srun}/portcullis/1-prep"
	log: SIM_DIR+"/{srun}/portcullis/portcullis.{srun}.2-junc.log"
	threads: THREADS
	message: "Using portcullis to extract junctions from: {input}"
	shell: "{params.load} {RUN_TIME} portcullis junc -o {params.outdir}/portcullis -t {threads} {params.prepdir} > {log} 2>&1"



rule tidy:
	input: 
		r1=rules.merge_fq1.output,
		r2=rules.merge_fq2.output,
		bam=rules.merge_bam.output,
		bed=rules.portcullis_junc.output.bed,
		tab=rules.portcullis_junc.output.tab
	output:
		bam=OUT_DIR+"/"+NAME+".{srun}.bam", 
		bed=OUT_DIR+"/"+NAME+".{srun}.bed",
		tab=OUT_DIR+"/"+NAME+".{srun}.tab",
		r1=OUT_DIR+"/"+NAME+".{srun}.R1.fq.gz",
		r2=OUT_DIR+"/"+NAME+".{srun}.R2.fq.gz"
	params: 
		bam="reads/sim/{srun}/spanki/"+NAME+".{srun}.bam", 
		bed="reads/sim/{srun}/portcullis/2-junc/portcullis.junctions.bed",
		tab="reads/sim/{srun}/portcullis/2-junc/portcullis.junctions.tab",
		r1="reads/sim/{srun}/spanki/"+NAME+".{srun}.R1.fq.gz",
		r2="reads/sim/{srun}/spanki/"+NAME+".{srun}.R2.fq.gz",
		chunkdir=SIM_DIR+"/{srun}/spanki/chunks"
	message: "Linking output files for run"
	shell: "ln -sf {params.bam} {output.bam} && ln -sf {params.bed} {output.bed} && ln -sf {params.tab} {output.tab} && ln -sf {params.r1} {output.r1} && ln -sf {params.r2} {output.r2} && touch -h {output.bam} {output.bed} {output.tab} {output.r1} {output.r2} && rm -rf {params.chunkdir}"
		
