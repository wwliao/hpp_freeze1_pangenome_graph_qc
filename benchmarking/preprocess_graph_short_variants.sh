#!/bin/bash
# Usage: preprocess_graph_short_variants.sh <VCF file>
VCF=$1
PREFIX=$(basename $VCF .vcf.gz)
CHROMS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"

bcftools view -e 'SVLEN>=50 | SVLEN<=-50 | SVTYPE="INV"' \
              -Oz -o ${PREFIX}.short.vcf.gz \
              ${VCF} \
    && bcftools index -t ${PREFIX}.short.vcf.gz \
    && bcftools view -r ${CHROMS} \
                     -Oz -o ${PREFIX}.short.chr1-22.vcf.gz \
                     ${PREFIX}.short.vcf.gz \
    && bcftools index -t ${PREFIX}.short.chr1-22.vcf.gz \
    && rm ${PREFIX}.short.vcf.gz*
