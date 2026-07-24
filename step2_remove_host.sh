#!/bin/bash

### Step2 QC: Remove host sequence

# Software
bowtie2=bowtie2
bowtie2_build=bowtie2-build
seqkit=seqkit

# Directories
CleanData=01_CleanData
RemoveHost=02_RemoveHost
RefDir=Ref

Ref=host_reference.fasta

Samples=$(cat sample.list)

mkdir -p ${RemoveHost} ${RefDir}

### Build bowtie2 index of host reference genome
${bowtie2_build} ${Ref} ${RefDir}/host_ref

for SampleID in ${Samples}
do
	### Remove host reads by bowtie2
	time -p ${bowtie2} \
		-I 200 -X 400 \
		--mm --end-to-end --sensitive \
		-x ${RefDir}/host_ref \
		-1 ${CleanData}/${SampleID}.clean.fq1.gz \
		-2 ${CleanData}/${SampleID}.clean.fq2.gz \
		-p 10 \
		--un-conc-gz ${RemoveHost}/${SampleID}.nohost.fq%.gz \
		--al-conc-gz ${RemoveHost}/${SampleID}.host.fq%.gz \
		-S ${RemoveHost}/${SampleID}.host.sam \
		-k 1 \
		2> ${RemoveHost}/${SampleID}.bowtie2.host.log

	### Count retained reads
	${seqkit} stats ${RemoveHost}/${SampleID}.nohost.fq1.gz ${RemoveHost}/${SampleID}.nohost.fq2.gz \
		>> ${RemoveHost}/nohost.fastq.stats.txt
done
