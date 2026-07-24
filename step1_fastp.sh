#!/bin/bash

### Step1 QC: RawData to CleanData

# Software
fastp=fastp

# Directories
RawData=00_RawData
CleanData=01_CleanData

Samples=$(cat sample.list)

mkdir -p ${RawData} ${CleanData}

for SampleID in ${Samples}
do
	cd ${CleanData}

	## QC: filter low-quality sequences by fastp v0.23.4
	time -p ${fastp} -w 10 \
		--length_required 120 \
		--detect_adapter_for_pe \
		-i ../${RawData}/${SampleID}_1.raw.fq.gz \
		-I ../${RawData}/${SampleID}_2.raw.fq.gz \
		-o ${SampleID}.clean.fq1.gz \
		-O ${SampleID}.clean.fq2.gz \
		-j ${SampleID}.fastp.json \
		-h ${SampleID}.fastp.html \
		-Y 30 -u 30 -q 20 -n 5 -y

	cd ..
done
