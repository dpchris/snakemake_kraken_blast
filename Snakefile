shell.executable("/bin/bash")
configfile: "config.yaml"


samples=["B2F05_frenchwastewater_spiked_17-D395_FR-META01_140517_NextSeq500PE80_TruSeq",
			"B2F06_italianair_spiked_BCW_IT-META01_050517_NextSeq500_Nextera",
			"B2F10_germanwastewater_spiked_17-D516_GE-META01_NextSeq500_TruSeq",
			"B2F11_swedishwastewater_spiked_17-D549_SE-META01_NextSeq500_TruSeq",
			"B2F12_germanwastewater_unspiked_17-D789_GE-META01_NextSeq500_TruSeq"]
			#"ATCC10987_R034-24","POU46","F_hispaniensis_R032-12","ANSES_56"]


#samples=["ATCC10987_R034-24","F_hispaniensis_R032-12","ANSES_56","POU46"]
strains=["anthracis","kurstaki","10987","hispaniensis"]
tools=["blast"]#kraken
strands=["R1","R2","single"]

rule all:
	input:
		expand("{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.sorted.bam.bai",tool=tools,strain=strains,sample=samples),
		#expand("{tool}_alignment_{strain}/{sample}_aln_control_paired_{strain}.bam",tool=["kraken","blast"],strain=strains,sample=samples)
		expand("fastqc/untrimmed/{sample}_{strand}_fastqc.html",strand=["R1","R2"],sample=samples),
		expand("fastqc/trimmed/{sample}_trimmed_{strand}_fastqc.html",strand=["R1","R2","single"],sample=samples),
		expand("krona_results/{sample}_trimmed_{strand}_cdb.html",sample=samples,strand=strands),
		#expand("blast_alignment_hispaniensis/{sample}_aln_merged_hispaniensis_consensus.fa",sample=samples)
		#expand("kraken_results/{sample}_trimmed_cdb_paired.txt",sample=samples)
		#expand("retained_reads/{tool}_alignment/{sample}_retained_{strand}_final_reads_{strain}_paired.fq.gz",tool=tools,sample=samples,strain=strains,strand=strands)
		#expand("megablast_results/{sample}_trimmed_{strain}_blast_output_{strand}_sorted.txt",strain=strains,strand=strands,sample=samples)

#index of ref fasta for depletion of ribosomal DNA (16S,23S)
rule bwa_index:
	input:
		"reference_genomes/{sample}_genome.fa"
	output:
		"reference_genomes/{sample}_genome.fa.amb",
		"reference_genomes/{sample}_genome.fa.ann",
		"reference_genomes/{sample}_genome.fa.bwt",
		"reference_genomes/{sample}_genome.fa.pac",
		"reference_genomes/{sample}_genome.fa.sa"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="bwa_index"
	shell :
		"bwa index {input}"

#fastqc of raw reads
rule fastqc:
	input:
		"/store/plateformes/I2BC/I2BC_B2FORENSICS/{sample}.fastq.gz"
	output:
		"fastqc/untrimmed/{sample}_fastqc.html",
		"fastqc/untrimmed/{sample}_fastqc.zip"
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="fastqc.{sample}"
	shell:
		"fastqc {input} -o fastqc/untrimmed/ -t {threads}"

#trimming of reads, bbduk advantage : remove adapters and is multithreaded
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
		"ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=28 threads=12 overwrite=true" #trim of bp with quality under 28

#compression of trimmed reads
rule compression_into_gzip:
	input:
		"fastq/{sample}.fq"
	output:
		protected("fastq/{sample}.fq.gz")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="compression_into_gzip.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

#fastqc of trimmed reads
rule fastqc2:
	input:
		"fastq/{sample}.fq.gz"
	output:
		"fastqc/trimmed/{sample}_fastqc.html",
		"fastqc/trimmed/{sample}_fastqc.zip"
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="fastqc2.{sample}"
	shell:
		"fastqc {input} -o fastqc/trimmed/ -t {threads}"

#kraken single analyse with custom database (optional)
rule kraken:
	input:
		kdb="/store/EQUIPES/LGBMB/dchristiany/kraken_databases/kraken_customdb_28-11-17_modif_kurstaki",
		reads="fastq/{sample}.fq.gz" #l'étape de rDNA depletion est devenu optionnelle avec le filtre samtools view q 25 (minimum de qualité de mapping)
	output:
		protected("kraken_results/{sample}_cdb.txt")
	conda:
		"envs/genomic.yaml"
	threads:20
	params:
		mem=360,
		jobname="kraken.{sample}"
	shell:
		"kraken --db {input.kdb} --fastq-input --gzip-compressed --threads {threads} --output {output} {input.reads}"

#kraken paired analyse with custom database
rule kraken_paired:
	input:
		kdb="/store/EQUIPES/LGBMB/dchristiany/kraken_databases/kraken_customdb_28-11-17_modif_kurstaki",
		fv="fastq/{sample}_R1.fq.gz",
		rv="fastq/{sample}_R2.fq.gz" 					
	output:
		protected("kraken_results/{sample}_cdb_paired.txt")
	conda:
		"envs/genomic.yaml"
	threads:20
	params:
		mem=360,
		jobname="kraken_paired.{sample}"
	shell:
		"kraken --db {input.kdb} --fastq-input --gzip-compressed --paired --threads {threads} --output {output} {input.fv} {input.rv}"

#krona representation of kraken results
rule krona:
	input:
		"kraken_results/{sample}.txt"
	output:
		"krona_results/{sample}.html"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="krona.{sample}"
	shell:
		"ktImportTaxonomy {input} -o {output} -q 2 -t 3 -k "	

#get reads_id identified as species of interest
rule extract_reads_id:
	input:
		tax="taxonomy_files/taxonomy_tree_{strain}.txt",
		kraken="kraken_results/{sample}_trimmed_{strand}_cdb.txt"
	output:
		"kraken_reads_id/{sample}_trimmed_{strain}/kraken_reads_id_{strand}.txt"
	threads:1
	params:
		mem=10,
		jobname="extract_reads_id.{sample}_trimmed_{strand}_{strain}"
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
		single="kraken_reads_id/{sample}_trimmed_{strain}/kraken_reads_id_{strand}.txt",
		paired="kraken_reads_id/{sample}_trimmed_{strain}/kraken_paired_reads_id.txt"
	output:
		"kraken_reads_id/{sample}_trimmed_{strain}/kraken_reads_id_single_{strand,R[0-9]}.txt"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=10,
		jobname="get_single_reads_id.{sample}_trimmed_{strain}.reads_id{strand}"
	shell:
		"cat {input.single} | grep -Fwf {input.paired} -v > {output} || true"

#get fastq for reads with two mates identified as species of interest
rule subseq_paired:
	input:
		reads="fastq/{sample}_trimmed_{strand}.fq.gz",
		reads_id="kraken_reads_id/{sample}_trimmed_{strain}/kraken_paired_reads_id.txt"
	output:
		"kraken_fastq/{sample}_trimmed_{strain}_{strand,(R1|R2)}.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_paired.{sample}_trimmed_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#get fastq for reads which lost their mate
rule subseq_new_single:
	input:
		reads="fastq/{sample}_trimmed_{strand}.fq.gz",
		reads_id="kraken_reads_id/{sample}_trimmed_{strain}/kraken_reads_id_single_{strand}.txt"
	output:
		temp("kraken_fastq/{sample}_trimmed_{strain}_single_{strand,R[0-9]}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_new_single.{sample}_trimmed_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#we merge all single reads identified as species of interest into a single fastq
rule merge_single_reads:
	input:
		expand("kraken_fastq/{{sample}}_trimmed_{{strain}}_{strand}.fq",strand=["single_R1","single_R2"])
	output:
		"kraken_fastq/{sample}_trimmed_{strain}_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=20,
		jobname="merge_single_reads.{sample}_{strain}"
	shell:
		"cat {input} > {output}"

#compression of fastq
rule compression_into_gzip3:
	input: 
		"kraken_fastq/{sample}.fq"
	output:
		protected("kraken_fastq/{sample}.fq.gz")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=60,
		jobname="compression_into_gzip3.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

#fastq conversion into fasta
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

#megablast on fasta (reads identified as species of interest by kraken)
rule megablast:
	input:
		"kraken_fasta/{sample}_{strand}.fa"
	output:
		protected("megablast_results/{sample}_blast_output_{strand}.txt")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="megablast.{sample}_{strand}"
	shell:
		"blastn -db nt -num_threads {threads} -max_target_seqs 1 -query {input} -out {output} -perc_identity 95 -qcov_hsp_perc 90 "
		"-outfmt '6 qseqid sacc stitle staxid score bitscore evalue pident nident mismatch qcovhsp'"

#filter blast output to keep best results with a least a bitscore of 50 and at most 2 SNP
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
			shell("cat {input} | grep thuringiensis | egrep 'kurstaki|Bc601|YWC2-8|YC-10|galleriae' |  sed 's/ /_/g' | awk '{{if ($10<2) print $0}}' | sort -k1,1 -u > {output} || true" )
		elif wildcards.strain == 'hispaniensis':
			shell("cat {input} | egrep 'hispaniensis|novicida 3523' |  sed 's/ /_/g' | awk '{{if ($10<2) print $0}}' | sort -k1,1 -u > {output} || true ")
		else :
			shell("cat {input} | grep {wildcards.strain} |  sed 's/ /_/g' | awk '{{if ($10<2) print $0}}' | sort -k1,1 -u -u > {output} || true" )

#extract reads id identified as species of interest from blast result
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

#get reads id for still paired end reads (both mates of paired end reads are identified as species of interest)
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

#get reads id of new single reads (only one reads of the two mates is identified as species of interest)
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
		"cat {input.single} | grep -Fwf {input.paired} -v > {output} || true"

#create fastq from paired reads_id
rule blast_subseq_paired:
	input:
		reads="kraken_fastq/{sample}_trimmed_{strain}_{strand}.fq.gz",
		reads_id="blast_reads_id/{sample}_trimmed_{strain}_blast_paired_reads_id.txt"
	output:
		"blast_fastq/{sample}_trimmed_{strain}_{strand,R[0-9]}.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="blast_subseq_paired.{sample}_trimmed_blast_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#create fastq from reads id of reads with lost mate
rule blast_subseq_new_single:
	input:
		reads="kraken_fastq/{sample}_trimmed_{strain}_{strand}.fq.gz",
		reads_id="blast_reads_id/{sample}_trimmed_{strain}_blast_sing_reads_id_{strand}.txt"
	output:
		temp("blast_fastq/{sample}_trimmed_sing_{strand,R[0-9]}_temp_{strain}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_new_single.{sample}_trimmed_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

rule blast_subseq_single:
	input:
		reads="kraken_fastq/{sample}_trimmed_{strain}_single.fq.gz",
		reads_id="blast_reads_id/{sample}_trimmed_{strain}_blast_reads_id_single.txt"
	output:
		temp("blast_fastq/{sample}_trimmed_single_temp_{strain}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_single.{sample}_trimmed_{strain}_single"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output}"

#concatenate into a single fastq new single reads from R1 or R2
rule blast_merge_single_reads:
	input:
		expand("blast_fastq/{{sample}}_trimmed_{strand}_temp_{{strain}}.fq",strand=['sing_R1','sing_R2','single'])
	output:
		"blast_fastq/{sample}_trimmed_{strain}_single.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=20,
		jobname="blast_merge_single_reads.{sample}_{strain}"
	shell:
		"cat {input} > {output}"

#compression of new fastqs
rule compression_into_gzip4:
	input: 
		"blast_fastq/{sample}.fq"
	output:
		protected("blast_fastq/{sample}.fq.gz")
	conda:
		"envs/genomic.yaml"
	threads:6
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
		index="reference_genomes/{strain}_genome.fa.amb",
		r1="{tool}_fastq/{sample}_trimmed_{strain}_R1.fq.gz",
		r2="{tool}_fastq/{sample}_trimmed_{strain}_R2.fq.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_paired_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		mem=60,
		jobname="bwa_map2.{tool}_{sample}_trimmed_{strain}"
	threads:6
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.r1} {input.r2} | "
		"samtools view -q 25 -@ {threads} -Sb -o {output})" #on ne garde pas les reads avec plus d'un mismatch | grep NM:i:[2-9] -v 


#mapping of single reads identified as species of interest
rule bwa_map_single:
	input:
		ref="reference_genomes/{strain}_genome.fa",
		index="reference_genomes/{strain}_genome.fa.amb",
		s="{tool}_fastq/{sample}_trimmed_{strain}_single.fq.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_single_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		jobname="bwa_map_single.{tool}_{sample}_trimmed_{strain}"
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.s} | "
		"samtools view -q 25 -@ {threads} -Sb -o {output})" #

#merge of the two previous alignment into a single
rule merge_bam:
	input:
		expand("{{tool}}_alignment_{{strain}}/{{sample}}_aln_{type}_{{strain}}.bam",type=["paired","single"])
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="merge_bam.{tool}_{sample}_{strain}"
	shell:
		"samtools merge -f -@ {threads} {output} {input}"

rule fixmate:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.bam"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.fixed.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="fixmate.{tool}_{sample}_{strain}"
	shell:
		"samtools fixmate -r {input} {output}"

rule sorting:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.fixed.bam"
	output:
		protected("{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.sorted.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="sorting.{tool}_{sample}_{strain}"
	shell:
		"samtools sort -@ threads {input} -o {output}"

rule index:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.sorted.bam"
	output:
		"{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.sorted.bam.bai"
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="index.{tool}_{sample}_{strain}"
	shell:
		"samtools index {input}"


#####alignement with reads with 1 mismatch max (control sample)

#mapping of paired-end reads (both mates identified as species of interest)
rule control_bwa_map2:
	input :
		ref="reference_genomes/{strain}_genome.fa",
		index="reference_genomes/{strain}_genome.fa.amb",
		r1="{tool}_fastq/{sample}_trimmed_{strain}_R1.fq.gz",
		r2="{tool}_fastq/{sample}_trimmed_{strain}_R2.fq.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_control_paired_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		mem=60,
		jobname="control_bwa_map2.{tool}_{sample}_trimmed_{strain}"
	threads:6
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.r1} {input.r2} | "
		"grep NM:i:[2-9] -v | samtools view -q 25 -@ {threads} -Sb -o {output})" #on ne garde pas les reads avec plus d'un mismatch | grep NM:i:[2-9] -v 


#mapping of single reads identified as species of interest
rule control_bwa_map_single:
	input:
		ref="reference_genomes/{strain}_genome.fa",
		index="reference_genomes/{strain}_genome.fa.amb",
		s="{tool}_fastq/{sample}_trimmed_{strain}_single.fq.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_control_single_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		rg="@RG\\tID:{sample}\\tSM:{sample}",
		jobname="control_bwa_map_single.{tool}_{sample}_trimmed_{strain}"
	shell:
		"(bwa mem -M -R '{params.rg}' -t {threads} {input.ref} {input.s} | "
		"grep NM:i:[2-9] -v | samtools view -q 25 -@ {threads} -Sb -o {output})" #on ne garde pas les reads avec plus d'un mismatch

#merge of the two previous alignment into a single
rule control_merge_bam:
	input:
		expand("{{tool}}_alignment_{{strain}}/{{sample}}_aln_control_{type}_{{strain}}.bam",type=["paired","single"])
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="control_merge_bam.{tool}_{sample}_{strain}"
	shell:
		"samtools merge -f -@ {threads} {output} {input}"

rule control_fixmate:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.bam"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.fixed.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="control_fixmate.{tool}_{sample}_{strain}"
	shell:
		"samtools fixmate -r {input} {output}"

rule control_sorting:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.fixed.bam"
	output:
		protected("{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.sorted.bam")
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="control_sorting.{tool}_{sample}_{strain}"
	shell:
		"samtools sort -@ threads {input} -o {output}"

rule control_index:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.sorted.bam"
	output:
		"{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.sorted.bam.bai"
	conda:
		"envs/genomic.yaml"
	threads:6
	params:
		mem=60,
		jobname="control_index.{tool}_{sample}_{strain}"
	shell:
		"samtools index {input}"


######## Variant call & consensus fasta


#identifying genomic variants with freebayes :
rule freebayes:
	input:
		bam="{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.sorted.bam",
		ref="reference_genomes/{strain}_genome.fa"
	output:
		temp("{tool}_alignment_{strain}/{sample}_aln_control_merged_{strain}.vcf")
	conda:
		"envs/variant_call.yaml"
	threads:1
	params:
		jobname="freebayes.{sample}_{strain}",
		mem=60
	shell:
		#freebayes -f <ref.fasta> -p 1 -i -u -X <aln.rmdup.sorted.bam> > var.vcf #-p 1 : haploide, -i : without indels, -u : without complex, -X without
		"freebayes -f {input.ref} -p 1 -i -u -X -m 25 {input.bam} | vcffilter -f 'QUAL > 100' > {output}"


#remove SNP if they are not at least at a distance of 12 bp from another SNP
#rule vcftools:
#	input:
#		"{tool}_alignment_{strain}/{sample}.vcf"
#	output:
#		"{tool}_alignment_{strain}/{sample}.filtered.vcf"
#	conda:
#		"envs/variant_call.yaml"
#	threads:1
#	params:
#		jobname="vcftools.{sample}_{strain}",
#		mem=60
#	shell:
#		"vcftools --vcf {input} --recode-INFO-all --thin 12 --out {output}"


rule bgzip:
	input:
		"{tool}_alignment_{strain}/{sample}.vcf"
	output:
		temp("{tool}_alignment_{strain}/{sample}.vcf.gz")
	conda:
		"envs/variant_call.yaml"
	threads:1
	params:
		jobname="bgzip.{sample}",
		mem=60
	shell:
		"bgzip -c {input} > {output}"

rule tabix:
	input:
		"{tool}_alignment_{strain}/{sample}.vcf.gz"
	output:
		temp("{tool}_alignment_{strain}/{sample}.vcf.gz.tbi")
	conda:
		"envs/variant_call.yaml"
	threads:1
	params:
		jobname="tabix.{sample}",
		mem=60
	shell:
		"tabix -p vcf {input}"

rule bcftools_consensus:
	input :
		ref="reference_genomes/{strain}_genome.fa",
		vcf="{tool}_alignment_{strain}/{sample}.vcf.gz",
		index="{tool}_alignment_{strain}/{sample}.vcf.gz.tbi"
	output :
		"{tool}_alignment_{strain}/{sample}_consensus.fa"
	conda:
		"envs/variant_call.yaml"
	threads:1
	params:
		jobname="bcftools_consensus.{sample}",
		mem=60
	shell:
		"bcftools consensus -f {input.ref} {input.vcf} > {output}"

##### extract retained reads (paired ends)


rule extract_final_reads_id:
	input:
		"{tool}_alignment_{strain}/{sample}_aln_merged_{strain}.sorted.bam"
	output:
		"retained_reads/{tool}_alignment/{sample}_final_reads_id_{strain}.txt"
	threads:1
	params:
		mem=10,
		jobname="extract_final_reads_id.{tool}_{sample}_{strain}"
	shell:
		"samtools view {input} | cut -f1 | grep \# -v > {output} || true"


rule final_reads_id_strand:
	input:
		reads="retained_reads/{tool}_alignment/{sample}_final_reads_id_{strain}.txt",
		reads_id="kraken_reads_id/{sample}_trimmed_{strain}/kraken_reads_id_R{strand}.txt"

	output:
		"retained_reads/{tool}_alignment/{sample}_final_reads_id_{strain}_R{strand}.txt"
	threads:1
	params:
		mem=10,
		jobname="final_reads_id_strand.{tool}_{sample}_{strain}"
	shell:
		"cat {input.reads} | grep -Fwf {input.reads_id} > {output} || true"

rule subseq_final_reads:
	input:
		reads="fastq/{sample}_trimmed_{strand}.fq.gz",
		reads_id="retained_reads/{tool}_alignment/{sample}_final_reads_id_{strain}_{strand}.txt"
	output:
		temp("retained_reads/{tool}_alignment/{sample}_retained_{strand}_final_reads_{strain}_{strand2,(R1|R2)}.fq")
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=60,
		jobname="subseq_final_reads.{sample}_{strain}_{strand}"
	shell:
		"seqtk subseq {input.reads} {input.reads_id} > {output} "


rule merge_fastq:
	input:
		reads1="retained_reads/{tool}_alignment/{sample}_retained_{strand}_final_reads_{strain}_R1.fq",
		reads2="retained_reads/{tool}_alignment/{sample}_retained_{strand}_final_reads_{strain}_R2.fq"
	output:
		"retained_reads/{tool}_alignment/{sample}_retained_{strand}_final_reads_{strain}_paired.fq"
	conda:
		"envs/genomic.yaml"
	threads:1
	params:
		mem=80,
		jobname="merge_fastq.{sample}_{strain}_{strand}"
	shell:
		"seqtk mergepe {input.reads1} {input.reads2} > {output} "

rule pigz_final_reads:
	input:
		"retained_reads/{tool}_alignment/{sample}.fq"
	output:
		protected("retained_reads/{tool}_alignment/{sample}.fq.gz")
	conda:
		"envs/genomic.yaml"
	threads:12
	params:
		mem=120,
		jobname="final_compression_into_gzip.{sample}"
	shell:
		"pigz -9 -p {threads} {input}"

