#!/bin/bash

### Step5 Genome binning and bin refinement

# Conda environments
MAGsEnv=MAGs
ConcoctEnv=concoct
ComebinEnv=comebin_env

# Software
bwa=bwa
samtools=samtools
jgi_summarize_bam_contig_depths=jgi_summarize_bam_contig_depths
metabat2=metabat2
cut_up_fasta=cut_up_fasta.py
concoct_coverage_table=concoct_coverage_table.py
concoct=concoct
merge_cutup_clustering=merge_cutup_clustering.py
run_comebin=run_comebin.sh

# Scripts
SplitConcoctBins=split_concoct_bins.py
BinRefinement=bin_refinement.sh

# Directories
RemoveHost=02_RemoveHost
Assembly=03_Megahit
Binning=05_Binning

Samples=$(cat sample.list)

mkdir -p ${Binning}

for SampleID in ${Samples}
do
	mkdir -p ${Binning}/${SampleID}
	cd ${Binning}/${SampleID}

	### Map clean reads back to assembled contigs
	source activate ${MAGsEnv}

	ln -sf ../../${Assembly}/${SampleID}/megahit/final.contigs.fa ${SampleID}.assembly.fa
	${bwa} index ${SampleID}.assembly.fa
	${bwa} mem -v 1 -t 10 ${SampleID}.assembly.fa \
		../../${RemoveHost}/${SampleID}.nohost.fq1.gz \
		../../${RemoveHost}/${SampleID}.nohost.fq2.gz \
		> ${SampleID}.sam

	${samtools} sort -@ 10 -O BAM -o ${SampleID}.sorted.bam ${SampleID}.sam
	${samtools} index -@ 10 -b ${SampleID}.sorted.bam
	rm -rf ${SampleID}.sam

	### Binning by MetaBAT2
	mkdir -p metabat2_bins
	cd metabat2_bins

	${jgi_summarize_bam_contig_depths} --outputDepth metabat_depth.txt ../${SampleID}.sorted.bam
	${metabat2} -i ../${SampleID}.assembly.fa \
		-a metabat_depth.txt \
		-o bin \
		-m 1500 \
		-t 10

	mkdir -p bin_result
	mv bin*.fa bin_result
	rm -rf bin_result/bin.lowDepth.fa bin_result/bin.tooShort.fa bin_result/bin.unbinned.fa
	cd ..

	source deactivate

	### Binning by CONCOCT
	source activate ${ConcoctEnv}
	mkdir -p concoct_bins
	cd concoct_bins

	${cut_up_fasta} ../${SampleID}.assembly.fa \
		-c 10000 \
		--merge_last \
		-b assembly_10K.bed \
		-o 0 \
		> assembly_10K.fa

	${concoct_coverage_table} assembly_10K.bed ../${SampleID}.sorted.bam > concoct_depth.txt
	${concoct} -l 1500 \
		-t 10 \
		--coverage_file concoct_depth.txt \
		--composition_file assembly_10K.fa \
		-b concoct_out

	${merge_cutup_clustering} concoct_out_clustering_gt1500.csv \
		> concoct_out_clustering_gt1500_merged.csv

	mkdir -p bin_result
	python ${SplitConcoctBins} concoct_out_clustering_gt1500_merged.csv ../${SampleID}.assembly.fa bin_result
	rm -rf bin_result/unbinned.fa
	cd ..

	source deactivate

	### Binning by ComeBin
	source activate ${ComebinEnv}
	mkdir -p comebin_bins
	cd comebin_bins

	${run_comebin} \
		-a ../${SampleID}.assembly.fa \
		-o comebin_result \
		-p ../ \
		-t 10

	mkdir -p bin_result
	cp -rf comebin_result/comebin_res/cluster_res/weight_seed_kmeans*_bins/*.fa bin_result
	cd ..

	source deactivate

	### Bin refinement by metaWRAP
	mkdir -p binning_refiner/input
	cd binning_refiner/input

	ln -sf ../../metabat2_bins/bin_result metabat2_bins
	ln -sf ../../concoct_bins/bin_result concoct_bins
	ln -sf ../../comebin_bins/bin_result comebin_bins
	cd ..

	${BinRefinement} \
		-o result.refinement \
		-A input/metabat2_bins \
		-B input/concoct_bins \
		-C input/comebin_bins \
		-t 10 \
		-m 70 \
		-c 50 \
		-x 10

	cd ../../..
done
