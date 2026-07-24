#!/bin/bash

### Step8 Salmon index and quantification

# Conda environment
SalmonEnv=salmon

# Software
salmon=salmon

# Scripts
SummarizeSalmon=lib/summarize_salmon_files.py
SplitSalmon=lib/split_salmon_out_into_bins.py

# Directories
ReferenceGenomes=dereplicated_genomes
Reads=02_RemoveHost
SalmonRoot=08_Salmon
RenameDir=rename

# Files and parameters
Reference=${SalmonRoot}/ref.fa
Index=${SalmonRoot}/ref.index
Samples=$(cat sample.list)
IndexThreads=20
QuantThreads=10

mkdir -p ${SalmonRoot}/result

### Build Salmon reference
cat ${ReferenceGenomes}/*.fa > ${Reference}

source activate ${SalmonEnv}
${salmon} index \
	-p ${IndexThreads} \
	-t ${Reference} \
	-i ${Index}

for SampleID in ${Samples}
do
	mkdir -p ${SalmonRoot}/result/${SampleID}
	cd ${SalmonRoot}/result/${SampleID}

	### Quantify paired-end reads by Salmon
	${salmon} quant \
		-i ../../ref.index \
		--libType IU \
		-1 ../../../${Reads}/${SampleID}.nohost.fq1.gz \
		-2 ../../../${Reads}/${SampleID}.nohost.fq2.gz \
		-o ${SampleID}.quant \
		--meta \
		-p ${QuantThreads}

	### Summarize and split Salmon results
	python ../../../${SummarizeSalmon}
	mkdir -p quant_files

	for f in *.quant.counts
	do
		mv ${f} quant_files/
	done

	python ../../../${SplitSalmon} \
		quant_files \
		../../../${RenameDir}/ \
		../../ref.fa \
		> bin_abundance_table.tab

	cd ../../..
done
