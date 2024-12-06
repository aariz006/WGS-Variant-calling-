# Germline Variant Calling Pipeline

This repository contains scripts for a germline variant calling pipeline using NGS data.

## Features
1. **Preprocessing**  
   - Download and prepare reference files.  
   - Perform quality control with `FastQC` and `MultiQC`.  
   - Align reads with `BWA` and process with `GATK`.  

2. **Variant Calling**  
   - Call variants with `GATK HaplotypeCaller`.  
   - Separate SNPs and indels into individual files.  

3. **Filtering and Annotation**  
   - Filter variants based on quality metrics.  
   - Annotate variants using `Funcotator`.  
   - Generate tabular summaries for downstream analysis.  

## Usage
```bash
bash wgs1.sh    # Run preprocessing
bash wgs2.sh    # Run filtering and annotation

## Dependencies
- **Tools**: `Samtools`, `BWA`, `GATK`, `FastQC`, `MultiQC`, `Funcotator`
- **Note**: Ensure proper paths are set in the scripts before running.



