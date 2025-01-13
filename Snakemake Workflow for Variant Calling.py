# Snakemake Workflow for Variant Calling with Explicit Software Paths
# shell behavior globally (optional)
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# working directory (optional, depending on your setup)
workdir: "/home/anum/cell_lines_analysis/"
rule all:
    input:
        "output/annotated_variants.vcf"

rule quality_control:
    input:
        fastq1="/home/anum/cell_lines_analysis/input/sample_1.fastq.gz",
        fastq2="/home/anum/cell_lines_analysis/input/sample_2.fastq.gz"
    output:
        "/home/anum/cell_lines_analysis/output/fastqc_done.txt"
    params:
        fastqc="/home/anum/software/fastqc/fastqc"  # Path to FastQC
    shell:
        "{params.fastqc} -o /home/anum/cell_lines_analysis/output {input.fastq1} {input.fastq2} && echo 'QC complete' > {output}"

rule trimming:
    input:
        fastq1="/home/anum/cell_lines_analysis/input/sample_1.fastq.gz",
        fastq2="/home/anum/cell_lines_analysis/input/sample_2.fastq.gz"
    output:
        trimmed1="/home/anum/cell_lines_analysis/output/trimmed_1_paired.fq.gz",
        trimmed2="/home/anum/cell_lines_analysis/output/trimmed_2_paired.fq.gz"
    params:
        trimmomatic="/home/anum/software/Trimmomatic-0.39/trimmomatic.jar",  # Path to Trimmomatic
        adapter_file="/home/anum/software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa"
    shell:
        "java -jar {params.trimmomatic} PE -phred33 {input.fastq1} {input.fastq2} \
         {output.trimmed1} /home/anum/cell_lines_analysis/output/trimmed_1_unpaired.fq.gz \
         {output.trimmed2} /home/anum/cell_lines_analysis/output/trimmed_2_unpaired.fq.gz \
         ILLUMINACLIP:{params.adapter_file}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"

rule alignment:
    input:
        fastq1="/home/anum/cell_lines_analysis/output/trimmed_1_paired.fq.gz",
        fastq2="/home/anum/cell_lines_analysis/output/trimmed_2_paired.fq.gz"
    output:
        bam="/home/anum/cell_lines_analysis/output/sorted.bam"
    params:
        bwa="/home/anum/software/bwa-0.7.17/bwa",  # Path to BWA
        samtools="/home/anum/software/samtools-1.15.1/samtools",  # Path to Samtools
        reference="/home/anum/software/datasets/GRCh38.fasta"
    shell:
        "{params.bwa} mem -t 8 {params.reference} {input.fastq1} {input.fastq2} | \
         {params.samtools} view -bS - | {params.samtools} sort -o {output.bam}"

rule mark_duplicates:
    input:
        bam="/home/anum/cell_lines_analysis/output/sorted.bam"
    output:
        bam="/home/anum/cell_lines_analysis/output/marked_duplicates.bam",
        metrics="/home/anum/cell_lines_analysis/output/marked_dup_metrics.txt"
    params:
        gatk="/home/anum/software/gatk-4.2/gatk"  # Path to GATK
    shell:
        "{params.gatk} MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} && \
         /home/anum/software/samtools-1.15.1/samtools index {output.bam}"

rule base_recalibration:
    input:
        bam="/home/anum/cell_lines_analysis/output/marked_duplicates.bam",
        known_sites="/home/anum/software/datasets/known_variants.vcf"
    output:
        recal_table="/home/anum/cell_lines_analysis/output/recal_data.table",
        recal_bam="/home/anum/cell_lines_analysis/output/recalibrated.bam"
    params:
        reference="/home/anum/software/datasets/GRCh38.fasta",
        gatk="/home/anum/software/gatk-4.2/gatk"  # Path to GATK
    shell:
        "{params.gatk} BaseRecalibrator -I {input.bam} -R {params.reference} \
         --known-sites {input.known_sites} -O {output.recal_table} && \
         {params.gatk} ApplyBQSR -R {params.reference} -I {input.bam} \
         --bqsr-recal-file {output.recal_table} -O {output.recal_bam}"

rule variant_calling:
    input:
        bam="/home/anum/cell_lines_analysis/output/recalibrated.bam"
    output:
        vcf="/home/anum/cell_lines_analysis/output/raw_variants.vcf"
    params:
        reference="/home/anum/software/datasets/GRCh38.fasta",
        gatk="/home/anum/software/gatk-4.2/gatk"  # Path to GATK
    shell:
        "{params.gatk} HaplotypeCaller -R {params.reference} -I {input.bam} -O {output.vcf}"

rule variant_filtering:
    input:
        vcf="/home/anum/cell_lines_analysis/output/raw_variants.vcf"
    output:
        vcf="/home/anum/cell_lines_analysis/output/filtered_variants.vcf"
    params:
        reference="/home/anum/software/datasets/GRCh38.fasta",
        gatk="/home/anum/software/gatk-4.2/gatk"  # Path to GATK
    shell:
        "{params.gatk} VariantFiltration -R {params.reference} -V {input.vcf} -O {output.vcf} \
         --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0' --filter-name 'basic_snp_filter'"

rule annotation:
    input:
        vcf="/home/anum/cell_lines_analysis/output/filtered_variants.vcf"
    output:
        vcf="/home/anum/cell_lines_analysis/output/annotated_variants.vcf"
    params:
        annovar="/home/anum/software/annovar/table_annovar.pl",  # Path to ANNOVAR
        annovar_db="/home/anum/software/annovar/humandb/",
        buildver="hg38"
    shell:
        "{params.annovar} {input.vcf} {params.annovar_db} -buildver {params.buildver} \
         -out /home/anum/cell_lines_analysis/output/annotated_variants \
         -remove -protocol refGene,cytoBand,dbsnp152,gnomad211_exome -operation g,r,f,f \
         -nastring . -vcfinput"
