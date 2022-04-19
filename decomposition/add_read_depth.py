#!/usr/bin/env python3
import re
import gzip
import argparse
from statistics import mean
from os.path import basename

import numpy as np
from cyvcf2 import VCF, Writer

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sample", required=True)
parser.add_argument("depth_matrix", help="matrix of read depth per edge")
parser.add_argument("vcffile", help="traversal-decomposed VCF file")
args = parser.parse_args()

def reverse_edge(edge):
    d = str.maketrans("><", "<>")
    step1, step2 = edge
    rev_edge = (step2.translate(d), step1.translate(d))
    return rev_edge

depth = {}
with gzip.open(args.depth_matrix, mode="rt") as infile:
    for i, line in enumerate(infile):
        cols = line.strip().split("\t")
        if i == 0:
            idx = cols.index(args.sample)
        else:
            if cols[1] == "+":
                step1 = f">{cols[0]}"
            else:
                step1 = f"<{cols[0]}"
            if cols[3] == "+":
                step2 = f">{cols[2]}"
            else:
                step2 = f"<{cols[2]}"
            depth[(step1, step2)] = int(cols[idx])

vcf = VCF(args.vcffile)

vcf.add_format_to_header({
    "ID": "DP",
    "Number": "1",
    "Type": "Integer",
    "Description": "Read depth"
})

vcf.add_format_to_header({
    "ID": "AD",
    "Number": "R",
    "Type": "Integer",
    "Description": "Read depth for each allele"
})

prefix = re.search("(\S+)\.vcf(?:\.gz)?$", basename(args.vcffile))[1]
w = Writer(f"{prefix}.added_depth.vcf.gz", vcf)
pattern = re.compile("([><]\d+)")
for v in vcf:
    traversals = v.INFO.get("AT").split(",")
    traversal_depths = []
    for traversal in traversals:
        steps = pattern.findall(traversal)
        depths = []
        for i in range(len(steps) - 1):
            edge = (steps[i], steps[i+1])
            if edge in depth:
                depths.append(depth[edge])
            else:
                depths.append(depth[reverse_edge(edge)])
        traversal_depths.append(int(mean(depths)))
    v.set_format("DP", np.array([sum(traversal_depths)]))
    v.set_format("AD", np.array(traversal_depths))
    w.write_record(v)
w.close()
vcf.close()
