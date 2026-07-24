#!/bin/bash

### Step3 Assembly: Assemble clean reads by MEGAHIT

# Software
megahit=megahit
seqkit=seqkit

# Directories
RemoveHost=02_RemoveHost
Assembly=03_Megahit

Samples=$(cat sample.list)

mkdir -p ${Assembly}

for SampleID in ${Samples}
do
	mkdir -p ${Assembly}/${SampleID}
	cd ${Assembly}/${SampleID}

	### Metagenomic assembly by MEGAHIT
	time -p ${megahit} \
		-1 ../../${RemoveHost}/${SampleID}.nohost.fq1.gz \
		-2 ../../${RemoveHost}/${SampleID}.nohost.fq2.gz \
		-o megahit \
		-t 10 \
		--min-contig-len 1500

	### Count assembled contigs
	${seqkit} stats megahit/final.contigs.fa \
		>> ../megahit.contig.stats.txt

	cd ../..
done
