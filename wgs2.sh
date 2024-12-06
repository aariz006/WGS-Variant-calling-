#!/bin/bash				

#for filtering and annotation	

# PATH AND DIRECTORIES

ref="/home/aariz/wgs/hg38/hg38.fa"
results="/home/aariz/wgs/results"


# FILTER VARIANTS

# STEP 1: FILTER SNPS

gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

# FILTER INDELS 

gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"


# SELECT THE VARIANT THAT PASS THE FILTER

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis_ready_snps.vcf 

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis_ready_indels.vcf 
	
	
	
# TO EXCLUDE VARIANT THAT FAILED THE FILTERS

cat analysis_ready_snps.vcf | grep -v -E "DP_filter|GQ_filter" > analysis_ready_snps_filteredGT.vcf
cat analysis_ready_indels.vcf | grep -v -E "DP_filter|GQ_filter" > analysis_ready_indels_filteredGT.vcf



# STEP 2 : VARIANT ANNOTATION(GAT4 FUNCOTATER)

gatk Funcotator \
	--variant ${results}/analysis_ready_snps_filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/aariz/funcotater/funcotator_dataSources.v1.8.hg38.20230908g \
	--output ${results}/analysis_ready_snps_filteredGT_funcotated.vcf \
	--output-file-format VCF 
	
gatk Funcotator \
	--variant ${results}/analysis_ready_indels_filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path /home/aariz/funcotater/funcotator_dataSources.v1.8.hg38.20230908g \
	--output ${results}/analysis_ready_indels_filteredGT_funcotated.vcf \
	--output-file-format VCF 

# EXTRACT FIELD FROM A VCF TO A TAB-DELIMATED TABLE

gatk VariantsToTable \
	-V  ${results}/analysis_ready_snps_filteredGT_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/output_snps.table
	
gatk VariantsToTable \
	-V  ${results}/analysis_ready_indels_filteredGT_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O ${results}/output_indels.table
	
# MANIPULATION OF FILES TO GENERATE A TSV FILE

cat analysis_ready_snps_filteredGT_funcotated.vcf | grep "Funcotation fields are:" | sed 's+|+\t+g' > output_curated_variant.txt
cat output_snps.table | cut -f5 | grep "NBPF1" | sed 's+|+\t+g' >> output_snps.table 

#can be carried out for both snps or indels


	
