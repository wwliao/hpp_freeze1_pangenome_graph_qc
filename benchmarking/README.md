# Benchmarking of graph variants (work in progress)

## Short variants

1. Preprocess truth calls generated from DeepVariant:

	```sh
	preprocess_deepvariant_vcf.sh $SAMPLE.GRCh38_no_alt.deepvariant.vcf.gz
	```

2. Preprocess traversal-decomposed variants in a pangenome graph:

	```sh
	preprocess_graph_short_variants.sh $SAMPLE.$GRAPH.decomposed.nonref.rmdup.vcf.gz
	```

3. Benchmark short variants

	```sh
	rtg vcfeval -b ${TRUTH} \
	            -c ${QUERY} \
	            -e ${CONF} \
	            -t ${SDF} \
	            -o ${DIR} \
	            -T ${THREADS} \
	            -m annotate \
	            --all-records \
	            --ref-overlap \
	            --no-roc
	```

4. Annotate `baseline.vcf.gz` and `calls.vcf.gz` generated from `rtg vcfeval`

	```sh
	PREFIX=$(basename $VCF .vcf.gz)
	echo -e '##INFO=<ID=Regions,Number=.,Type=String,Description="Tags for regions">' > hdr.txt
	bcftools annotate -a stratification_regions.bed.gz \
	                  -c CHROM,FROM,TO,INFO/Regions \
	                  -h hdr.txt \
	                  -l Regions:unique \
	                  -Oz -o ${PREFIX}.annotated.vcf.gz \
	                  ${VCF} \
	    && bcftools index -t ${PREFIX}.annotated.vcf.gz
	```

5. Calculate performance metrics

	```sh
	quantify.py -o $SAMPLE.$GRAPH.performance.txt \
	            baseline.annotated.vcf.gz \
	            calls.annotated.vcf.gz
	```

## Structural variants

## Thoughts

1. Do we need to add `--Xloose-match-distance` (seems this option is only for GA4GH mode) when running `rtg vcfeval`?
2. Should we move the dipcall confident regions to the annotation step instead of feeding into `rtg vcfeval` at the benchmarking step?
