#!/bin/bash

### Step4 Taxonomic and functional profiling: MetaPhlAn and HUMAnN

# Conda environment
CondaEnv=biobakery

# Software
metaphlan=metaphlan
humann=humann
humann_join_tables=humann_join_tables
humann_renorm_table=humann_renorm_table
humann_regroup_table=humann_regroup_table
humann_rename_table=humann_rename_table

# Directories
RemoveHost=02_RemoveHost
MetaPhlAn=04_MetaPhlAn
HUMAnN=05_HUMAnN

Samples=$(cat sample.list)

mkdir -p ${MetaPhlAn} ${HUMAnN}

### Activate biobakery environment
source activate ${CondaEnv}

for SampleID in ${Samples}
do
	mkdir -p ${MetaPhlAn}/${SampleID}
	cd ${MetaPhlAn}/${SampleID}

	### Taxonomic profiling by MetaPhlAn
	time -p ${metaphlan} \
		../../${RemoveHost}/${SampleID}.nohost.fq1.gz,../../${RemoveHost}/${SampleID}.nohost.fq2.gz \
		--bowtie2out ${SampleID}.metagenome.bowtie2.bz2 \
		--nproc 10 \
		--input_type fastq \
		--unclassified_estimation \
		-o ${SampleID}.profiled_metagenome.txt

	cd ../..
done

for SampleID in ${Samples}
do
	mkdir -p ${HUMAnN}/${SampleID}
	cd ${HUMAnN}/${SampleID}

	### Functional profiling by HUMAnN
	time -p ${humann} \
		--threads 10 \
		--input ../../${RemoveHost}/${SampleID}.nohost.fq1.gz \
		--output .

	GeneFamilies=$(find . -maxdepth 1 -type f -name '*genefamilies.tsv' | head -n 1)
	${humann_renorm_table} \
		--input ${GeneFamilies} \
		--output ../${SampleID}_genefamilies-cpm.tsv \
		--units cpm \
		--update-snames

	cd ../..
done

### Join, regroup, renorm and rename HUMAnN tables
${humann_join_tables} \
	--input ${HUMAnN} \
	--output ${HUMAnN}/comb_genefamilies-cpm.tsv \
	--file_name genefamilies-cpm.tsv

${humann_regroup_table} \
	--input ${HUMAnN}/comb_genefamilies-cpm.tsv \
	--output ${HUMAnN}/comb_genefamilies-cpm.uniref90_ko.tsv \
	--groups uniref90_ko

${humann_renorm_table} \
	--input ${HUMAnN}/comb_genefamilies-cpm.uniref90_ko.tsv \
	--output ${HUMAnN}/comb_genefamilies-cpm.uniref90_ko.normed.tsv \
	--units cpm

${humann_rename_table} \
	--input ${HUMAnN}/comb_genefamilies-cpm.uniref90_ko.normed.tsv \
	--output ${HUMAnN}/comb_genefamilies-cpm.uniref90_ko.normed.named.tsv \
	--names kegg-orthology
