#!/usr/bin/env python3
import argparse
from collections import Counter
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--superpop", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-m", "--method", required=True)
parser.add_argument("vcffile")
args = parser.parse_args()

chroms = set(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"])

vcf = VCF(args.vcffile)
counter = Counter()
for v in vcf:
    if v.CHROM in chroms:
        vtype = v.INFO.get("SVTYPE")
        if vtype in ["SNP", "MNP"]:
            vtype = "SNP/MNP"
        if vtype in ["INS", "DEL"]:
            vlen = v.INFO.get("SVLEN")
            if -50 < vlen < 50:
                vtype = "INDEL"
        counter.update([vtype])
vcf.close()

with open(f"{args.sample}.{args.method}.vcount.txt", "w") as outfile:
    outfile.write("Superpopulation\tSample\tMethod\tVariantType\tCount\n")
    for vtype, count in counter.most_common():
        outfile.write(f"{args.superpop}\t{args.sample}\t{args.method}\t{vtype}\t{count}\n")
