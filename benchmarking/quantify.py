#!/usr/bin/env python3
import argparse
from collections import defaultdict
from cyvcf2 import VCF

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--outfile", required=True)
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("-g", "--graph", required=True)
parser.add_argument("truthvcf")
parser.add_argument("queryvcf")
args = parser.parse_args()

truth_count = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
vcf = VCF(args.truthvcf)
for v in vcf:
    decision = v.INFO.get("BASE")
    if decision in ["TP", "FN", "FN_CA"]:
        if len(v.REF) == len(v.ALT[0]):
            vtype = "SNP"
        else:
            vtype = "INDEL"
        truth_count["Confident"][vtype][decision] += 1
        truth_count["Confident"]["Both"][decision] += 1
        if v.INFO.get("Regions"):
            for region in v.INFO.get("Regions").split(","):
                truth_count[region][vtype][decision] += 1
                truth_count[region]["Both"][decision] += 1
vcf.close()

query_count = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
vcf = VCF(args.queryvcf)
for v in vcf:
    decision = v.INFO.get("CALL")
    if decision in ["TP", "FP", "FP_CA"]:
        if v.INFO.get("SVTYPE") in ["SNP", "MNP"]:
            vtype = "SNP"
        else:
            vtype = "INDEL"
        query_count["Confident"][vtype][decision] += 1
        query_count["Confident"]["Both"][decision] += 1
        if v.INFO.get("Regions"):
            for region in v.INFO.get("Regions").split(","):
                query_count[region][vtype][decision] += 1
                query_count[region]["Both"][decision] += 1
vcf.close()

regions = sorted(set(query_count) | set(truth_count))
vtypes = ["Both", "SNP", "INDEL"]

# Calculate recall, precision, and F1 score
with open(args.outfile, "w") as outfile:
    outfile.write("Sample\tGraph\tRegion\tVariant_Type\tTP-base\tTP-call\t")
    outfile.write("FN\tFN_CA\tFP\tFP_CA\tRecall\tPrecision\tF1\n")
    for region in regions:
        for vtype in vtypes:
            truth_tp = truth_count[region][vtype]["TP"]
            truth_fn = truth_count[region][vtype]["FN"]
            truth_fn_ca = truth_count[region][vtype]["FN_CA"]
            recall = truth_tp / (truth_tp + truth_fn + truth_fn_ca) * 100
            query_tp = query_count[region][vtype]["TP"]
            query_fp = query_count[region][vtype]["FP"]
            query_fp_ca = query_count[region][vtype]["FP_CA"]
            precision = query_tp / (query_tp + query_fp + query_fp_ca) * 100
            f1 = 2 * recall * precision / (recall + precision)
            outfile.write(f"{args.sample}\t{args.graph}\t{region}\t{vtype}\t{truth_tp}\t{query_tp}\t")
            outfile.write(f"{truth_fn}\t{truth_fn_ca}\t")
            outfile.write(f"{query_fp}\t{query_fp_ca}\t")
            outfile.write(f"{recall:.2f}\t{precision:.2f}\t{f1:.2f}\n")
