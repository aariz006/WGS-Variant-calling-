#!/bin/bash

# DOWNLOAD DATA

wget -P /home/aariz/wgs/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P /home/aariz/wgs/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

# DOWNLOAD REFERENCE FILE 

wget -P /home/aariz/wgs/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /home/aariz/wgs/hg38/hg38.fa

# index ref - .fia file before running haplotype caller 
samtools faidx /home/aariz/wgs/hg38/hg38.fa

# ref - dictionary file 

gatk CreateSequenceDictionary R=/home/aariz/wgs/hg38/hg38.fa O=/home/aariz/wgs/hg38/hg38.dict

# download known sites files for BSQR from GATK resource bundle

wget -P /home/aariz/wgs/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /home/aariz/wgs/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# path and directories

ref="/home/aariz/wgs/hg38/hg38.fa"
known_sites="/home/aariz/wgs/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/aariz/wgs/aligned_reads"
reads="/home/aariz/wgs/reads"
results="/home/aariz/wgs/results"
data="/home/aariz/wgs/data"
# VA
# Step 1 : QC 

echo "Step 1 : QC"
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/
multiqc ${reads}/*fastqc.zip -o ${reads}

# Step 2 : BWA INDEX REFERENCE

echo "BWA index reference"
bwa index ${ref}

# Step 3 : BWA ALIGNMENT

echo "Alighment using bwa"
bwa mem -t 5 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

# Step 4 : MAKR DUPLICATE READS AND SORE - gatk4

echo "Flag and sort reads"
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

# Step 5 : BASE QUALITY RECALIBRATION

# 1. build the model
echo "build the model"
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
echo "apply the modal"
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 

# Step 6 : QC-collect alignment and insert size metrices

echo "summarising qc"
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


# Step 7 : VARIANT CALLING - gatk haplotype caller 

echo "variant calling"
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

# Step 8 : EXTRACT SNPS and INDELS

echo "extracting snps and indels"
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
