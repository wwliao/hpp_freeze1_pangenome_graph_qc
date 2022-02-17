#!/bin/bash
# Usage: preprocess_deepvariant_vcf.sh <VCF file>
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

bcftools norm -m -any -Ou ${VCF} \
	| bcftools view -e 'GT="ref" | GT="0/." | GT="./0" | GT="./." | STRLEN(ALT)-STRLEN(REF)>=50 | STRLEN(ALT)-STRLEN(REF)<=-50' \
	                -f 'PASS' \
					-Oz -o ${PREFIX}.short.passed.vcf.gz \
    && bcftools index -t ${PREFIX}.short.passed.vcf.gz \
    && bcftools view -r ${CHROMS} \
	                 -Oz -o ${PREFIX}.short.passed.chr1-22.vcf.gz \
                     ${PREFIX}.short.passed.vcf.gz \
    && bcftools index -t ${PREFIX}.short.passed.chr1-22.vcf.gz \
    && rm ${PREFIX}.short.passed.vcf.gz*
