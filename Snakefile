#shell.prefix("source /home/david.christiany/soft/miniconda3/bin/activate genomic ;")
shell.executable("/bin/bash")
configfile: "config.yaml"


#samples=["B2F05_frenchwastewater_spiked_17-D395_FR-META01_140517_NextSeq500PE80_TruSeq",
#			"B2F06_italianair_spiked_BCW_IT-META01_050517_NextSeq500_Nextera",
#			"B2F10_germanwastewater_spiked_17-D516_GE-META01_NextSeq500_TruSeq",
#			"B2F11_swedishwastewater_spiked_17-D549_SE-META01_NextSeq500_TruSeq",
#			"B2F12_germanwastewater_unspiked_17-D789_GE-META01_NextSeq500_TruSeq"]

samples=["B2F12_germanwastewater_unspiked_17-D789_GE-META01_NextSeq500_TruSeq"]
strains=["anthracis","kurstaki","10987"]

#samples=config["samples"]

rule all:
	input:
		expand("{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.sorted.bam.bai",tool=["kraken","blast"],strain=strains,sample=samples),
		#expand("fastqc/untrimmed/{sample}_{strand}_fastqc.html",strand=["R1","R2"],sample=samples),
		expand("fastqc/trimmed/{sample}_trimmed_{strand}_fastqc.html",strand=["R1","R2","single"],sample=samples),
		expand("krona_results/{sample}_rDNA_depleted_{strand}_cdb.html",strand=["R1","R2","single"],sample=samples)

#index of ref fasta for depletion of ribosomal DNA (16S,23S)
rule bwa_index:
	input:
		"silva/silva_128_lsu_ssu_tax.fasta"
	output:
		"silva/silva_128_lsu_ssu_tax.fasta.amb",
		"silva/silva_128_lsu_ssu_tax.fasta.ann",
		"silva/silva_128_lsu_ssu_tax.fasta.bwt",
		"silva/silva_128_lsu_ssu_tax.fasta.pac",
		"silva/silva_128_lsu_ssu_tax.fasta.sa"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="bwa_index"
	shell :
		"bwa index {input}"

rule fastqc:
	input:
		"/store/plateformes/I2BC/I2BC_B2FORENSICS/{sample}.fastq.gz"
	output:
		"fastqc/untrimmed/{sample}_fastqc.html",
		"fastqc/untrimmed/{sample}_fastqc.zip"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="fastqc.{sample}"
	shell:
		"fastqc {input} -o fastqc/untrimmed/ -t {threads}"

rule bbduk:
	input:
		r1="/store/plateformes/I2BC/I2BC_B2FORENSICS/{sample}_R1.fastq.gz",
		r2="/store/plateformes/I2BC/I2BC_B2FORENSICS/{sample}_R2.fastq.gz"
	output:
		r1="fastq/{sample}_trimmed_R1.fq",
		r2="fastq/{sample}_trimmed_R2.fq",
		s="fastq/{sample}_trimmed_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="bbduk.{sample}"
	shell:
		"bbduk.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.s} "
		"ref=adapters.fa " #trim of adapters
		"ktrim=r k=23 mink=11 hdist=1 tpe tbo ftl=10 qtrim=r trimq=28 threads=12 overwrite=true" #trim of 10 first bp and quality under 28

rule compression_into_gzip:
	input:
		"fastq/{sample}.fq"
	output:
		"fastq/{sample}.fq.gz"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="compression_into_gzip.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

rule fastqc2:
	input:
		"fastq/{sample}.fq.gz"
	output:
		"fastqc/trimmed/{sample}_fastqc.html",
		"fastqc/trimmed/{sample}_fastqc.zip"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="fastqc2.{sample}"
	shell:
		"fastqc {input} -o fastqc/trimmed/ -t {threads}"

#mapping for rDNA_depletion
rule bwa_map:
	input :
		ref="silva/silva_128_lsu_ssu_tax.fasta",
		r1="fastq/{sample}_trimmed_R1.fq.gz",
		r2="fastq/{sample}_trimmed_R2.fq.gz"
	output:
		temp("rDNA_depleted/{sample}.bam")
	conda:
		"envs/genomic.yaml"
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		mem=120,
		jobname="bwa_map.{sample}"
	threads: 12
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.r1} {input.r2} | "
		"samtools view -@ {threads} -Sb -o {output})"

rule samtools_sort:
	input:
		"rDNA_depleted/{sample}.bam"
	output:
		"rDNA_depleted/{sample}_sorted.bam"
	conda:
		"envs/genomic.yaml"
	threads: 12
	params:
		mem=60,
		jobname="samtools_sort.{sample}"
	shell:
		"samtools sort -n -@ {threads} {input} -o {output}"

#extraction of unmapped reads on ribosomal DNA
rule reads_extraction:
	input:
		"rDNA_depleted/{sample}_sorted.bam"
	output:
		ffastq="rDNA_depleted/{sample}_rDNA_depleted_R1.fq",
		rfastq="rDNA_depleted/{sample}_rDNA_depleted_R2.fq",
		sfastq="rDNA_depleted/{sample}_rDNA_depleted_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="reads_extraction.{sample}"
	shell:
		"samtools bam2fq -@ {threads} -f4 {input} -1 {output.ffastq} -2 {output.rfastq} -s {output.sfastq} -n "

rule compression_into_gzip2:
	input: 
		"rDNA_depleted/{sample}.fq"
	output:
		"rDNA_depleted/{sample}.fq.gz"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=60,
		jobname="compression_into_gzip2.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

rule kraken:
	input:
		kdb="/store/EQUIPES/LGBMB/dchristiany/kraken_databases/kraken_customdb_28-11-17_modif_kurstaki",
		reads="rDNA_depleted/{sample}.fq.gz"
	output:
		"kraken_results/{sample}_cdb.txt"
	conda:
		"envs/genomic.yaml"
	threads:20
	params:
		mem=360,
		jobname="kraken.{sample}"
	shell:
		"kraken --db {input.kdb} --fastq-input --gzip-compressed --threads {threads} --output {output} {input.reads}"

rule krona:
	input:
		"kraken_results/{sample}.txt"
	output:
		"krona_results/{sample}.html"
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="krona.{sample}"
	shell:
		"ktImportTaxonomy {input} -o {output} -q 2 -t 3 -s 0"	

#get reads_id identified as species of interest
rule extract_reads_id:
	input:
		tax="taxonomy_files/taxonomy_tree_{strain}.txt",
		kraken="kraken_results/{sample}_rDNA_depleted_{strand}_cdb.txt"
	output:
		"kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_reads_id_{strand}.txt"
	threads:1
	params:
		mem=10,
		jobname="extract_reads_id.{sample}_rDNA_depleted_{strand}_{strain}"
	shell:
		"cat {input.kraken} | grep -Ff {input.tax} | cut -f2 | sort > {output} || true" #||true avoid non zero status exit 1 error with empty result file

#get reads_id for reads with two mates identified as species of interest
rule get_paired_reads_id:
	input:
		expand("kraken_reads_id/{{sample}}/kraken_reads_id_{strand}.txt",strand=["R1","R2"])
	output:
		"kraken_reads_id/{sample}/kraken_paired_reads_id.txt"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="get_paired_reads_id.{sample}"
	shell:
		"comm -12 {input} > {output}"

#get reads_id for reads which lost their mate
rule get_single_reads_id:
	input:
		single="kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_reads_id_{strand}.txt",
		paired="kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_paired_reads_id.txt"
	output:
		"kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_reads_id_single_{strand,R[0-9]}.txt"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="get_single_reads_id.{sample}_rDNA_depleted_{strain}.reads_id{strand}"
	shell:
		"cat {input.single} | grep -Ff {input.paired} -v > {output} || true"

#get fastq for reads with two mates identified as species of interest
rule subseq_paired:
	input:
		reads="rDNA_depleted/{sample}_rDNA_depleted_{strand}.fq.gz",
		reads_id="kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_paired_reads_id.txt"
	output:
		"kraken_fastq/{sample}_rDNA_depleted_{strain}_{strand,(R1|R2)}.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_paired.{sample}_rDNA_depleted_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#get fastq for reads which lost their mate
rule subseq_new_single:
	input:
		reads="rDNA_depleted/{sample}_rDNA_depleted_{strand}.fq.gz",
		reads_id="kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_reads_id_single_{strand}.txt"
	output:
		temp("kraken_fastq/{sample}_rDNA_depleted_{strain}_single_{strand,R[0-9]}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_new_single.{sample}_rDNA_depleted_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#get fastq for single reads identified as species of interest
rule subseq_single:
	input:
		reads="rDNA_depleted/{sample}_rDNA_depleted_single.fq.gz",
		reads_id="kraken_reads_id/{sample}_rDNA_depleted_{strain}/kraken_reads_id_single.txt"
	output:
		temp("kraken_fastq/{sample}_rDNA_depleted_{strain}_temp_single.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_single.{sample}_rDNA_depleted_{strain}_single"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#we merge all single reads identified as species of interest into a single fastq
rule merge_single_reads:
	input:
		expand("kraken_fastq/{{sample}}_rDNA_depleted_{{strain}}_{strand}.fq",strand=["single_R1","single_R2","temp_single"])
	output:
		"kraken_fastq/{sample}_rDNA_depleted_{strain}_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=20,
		jobname="merge_single_reads.{sample}_{strain}"
	shell:
		"cat {input} > {output}"

rule compression_into_gzip3:
	input: 
		"kraken_fastq/{sample}.fq"
	output:
		"kraken_fastq/{sample}.fq.gz"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=60,
		jobname="compression_into_gzip3.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

rule fq_to_fa:
	input:
		"kraken_fastq/{sample}.fq.gz"
	output:
		"kraken_fasta/{sample}.fa"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="fq_to_fa.{sample}"
	shell:
		"seqtk seq -A {input} > {output}"

rule megablast:
	input:
		"kraken_fasta/{sample}_{strand}.fa"
	output:
		"megablast_results/{sample}_blast_output_{strand}.txt"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="megablast.{sample}_{strand}"
	shell:
		"blastn -db nt -num_threads {threads} -max_target_seqs 1 -query {input} -out {output} "
		"-outfmt '6 qseqid sacc stitle staxid score bitscore evalue'"

rule sort_blast_output:
	input:
		"megablast_results/{sample}_{strain}_blast_output_{strand}.txt"
	output:
		"megablast_results/{sample}_{strain}_blast_output_{strand}_sorted.txt"
	threads:1
	params:
		mem=10,
		jobname="sort_blast_output.{sample}_{strain}_{strand}"
	run:
		if wildcards.strain == 'kurstaki' :
			shell("cat {input} | grep thuringiensis | egrep 'kurstaki|Bc601|YWC2-8|YC-10' | sort -k1,1 -u > {output} || true" )
		else :
			shell("cat {input} grep {wildcards.strain} | sort -k1,1 -u > {output} || true" )

rule blast_reads_id: 
	input:
		"megablast_results/{sample}_blast_output_{strand}_sorted.txt"
	output:
		"blast_reads_id/{sample}_blast_reads_id_{strand}.txt"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="blast_reads_id.{sample}_{strand}"
	shell:
		"cat {input} | cut -f1 > {output}"

rule blast_reads_id_paired:
	input:
		expand("blast_reads_id/{{sample}}_blast_reads_id_{strand}.txt",strand=["R1","R2"])
	output:
		"blast_reads_id/{sample}_blast_paired_reads_id.txt"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="blast_reads_id_paired.{sample}"
	shell:
		"comm -12 {input} > {output}"


rule blast_reads_id_single:
	input:
		single="blast_reads_id/{sample}_blast_reads_id_{strand}.txt",
		paired="blast_reads_id/{sample}_blast_paired_reads_id.txt"
	output:
		"blast_reads_id/{sample}_blast_sing_reads_id_{strand}.txt" 
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="blast_reads_id_single.{sample}.reads_id_{strand}"
	shell:
		"cat {input.single} | grep -Ff {input.paired} -v > {output} || true"


rule blast_subseq_paired:
	input:
		reads="kraken_fastq/{sample}_rDNA_depleted_{strain}_{strand}.fq.gz",
		reads_id="blast_reads_id/{sample}_rDNA_depleted_{strain}_blast_reads_id_{strand}.txt"
	output:
		"blast_fastq/{sample}_rDNA_depleted_{strain}_{strand,R[0-9]}.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="blast_subseq_paired.{sample}_rDNA_depleted_blast_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

rule blast_subseq_new_single:
	input:
		reads="kraken_fastq/{sample}_rDNA_depleted_{strain}_{strand}.fq.gz",
		reads_id="blast_reads_id/{sample}_rDNA_depleted_{strain}_blast_sing_reads_id_{strand}.txt"
	output:
		temp("blast_fastq/{sample}_rDNA_depleted_sing_{strand,R[0-9]}_temp_{strain}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_new_single.{sample}_rDNA_depleted_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

rule blast_subseq_single:
	input:
		reads="kraken_fastq/{sample}_rDNA_depleted_{strain}_single.fq.gz",
		reads_id="blast_reads_id/{sample}_rDNA_depleted_{strain}_blast_reads_id_single.txt"
	output:
		temp("blast_fastq/{sample}_rDNA_depleted_single_temp_{strain}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_single.{sample}_rDNA_depleted_{strain}_single"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

rule blast_merge_single_reads:
	input:
		expand("blast_fastq/{{sample}}_rDNA_depleted_{strand}_temp_{{strain}}.fq",strand=['sing_R1','sing_R2','single'])
	output:
		"blast_fastq/{sample}_rDNA_depleted_{strain}_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=20,
		jobname="blast_merge_single_reads.{sample}_{strain}"
	shell:
		"cat {input} > {output}"

rule compression_into_gzip4:
	input: 
		"blast_fastq/{sample}.fq"
	output:
		"blast_fastq/{sample}.fq.gz"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=60,
		jobname="compression_into_gzip4.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

###alignments 

#mapping of paired-end reads (both mates identified as species of interest)
rule bwa_map2:
	input :
		ref="reference_genomes/{strain}_genome.fa",
		r1="{tool}_fastq/{sample}_rDNA_depleted_{strain}_R1.fq.gz",
		r2="{tool}_fastq/{sample}_rDNA_depleted_{strain}_R2.fq.gz",

	output:
		temp("{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_paired_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		mem=120,
		jobname="bwa_map2.{tool}_{sample}_rDNA_depleted_{strain}"
	threads: 12
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.r1} {input.r2} | "
		"samtools view -@ {threads} -Sb -o {output})"

#mapping of single reads identified as species of interest
rule bwa_map_single:
	input:
		ref="reference_genomes/{strain}_genome.fa",
		s="{tool}_fastq/{sample}_rDNA_depleted_{strain}_single.fq.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_single_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		jobname="bwa_map_single.{tool}_{sample}_rDNA_depleted_{strain}"
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.s} | "
		"samtools view -@ {threads} -Sb -o {output})"

#merge of the two previous alignment into a single
rule merge_bam:
	input:
		expand("{{tool}}_alignment_{{strain}}/{{sample}}_rDNA_depleted_{{strain}}_aln/{{sample}}_aln_{type}_{{strain}}.bam",type=["paired","single"])
	output:
		temp("{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="merge_bam.{tool}_{sample}_{strain}"
	shell:
		"samtools merge -f -@ {threads} {output} {input}"

rule fixmate:
	input:
		"{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.bam"
	output:
		temp("{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.fixed.bam")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="fixmate.{tool}_{sample}_{strain}"
	shell:
		"samtools fixmate -r {input} {output}"

rule sorting:
	input:
		"{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.fixed.bam"
	output:
		"{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.sorted.bam"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="sorting.{tool}_{sample}_{strain}"
	shell:
		"samtools sort -@ threads {input} -o {output}"

rule index:
	input:
		"{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.sorted.bam"
	output:
		"{tool}_alignment_{strain}/{sample}_rDNA_depleted_{strain}_aln/{sample}_aln_merged_{strain}.sorted.bam.bai"
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="index.{tool}_{sample}_{strain}"
	shell:
		"samtools index {input}"
