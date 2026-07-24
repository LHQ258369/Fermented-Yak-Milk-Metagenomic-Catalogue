#!/bin/bash

### Step9 Biosynthetic gene cluster prediction and clustering

# Conda environments
AntiSMASHEnv=antismash
BigSCAPEEnv=bigscape

# Software
antismash=antismash
bigscape=bigscape

# Files and directories
MAGList=MAG.list
AntiSMASHRoot=09_AntiSMASH
AntiSMASHOut=${AntiSMASHRoot}/antismash_out
BigSCAPEOut=${AntiSMASHRoot}/bigscape2_out_mibig4
PfamDB=pfam/Pfam-A.hmm

# Parameters
ParallelJobs=10
AntiSMASHCpus=3
BigSCAPECpus=20

mkdir -p ${AntiSMASHOut}/logs

run_antismash_one_mag() {
	Mag=$1
	Name=$(basename ${Mag})
	Name=${Name%.*}
	MagOut=${AntiSMASHOut}/${Name}
	LogFile=${AntiSMASHOut}/logs/${Name}.log

	conda run -n ${AntiSMASHEnv} ${antismash} \
		--taxon bacteria \
		--genefinding-tool prodigal \
		--cb-knownclusters \
		--cc-mibig \
		--fullhmmer \
		--cpus ${AntiSMASHCpus} \
		--output-dir ${MagOut} \
		--output-basename ${Name} \
		--logfile ${LogFile} \
		-v \
		${Mag}
}

export AntiSMASHEnv antismash AntiSMASHOut AntiSMASHCpus
export -f run_antismash_one_mag

### Biosynthetic gene cluster prediction by antiSMASH
grep -v '^[[:space:]]*$' ${MAGList} | \
	grep -v '^[[:space:]]*#' | \
	xargs -r -d '\n' -n 1 -P ${ParallelJobs} bash -c 'run_antismash_one_mag "$1"' _

### Biosynthetic gene cluster clustering by BiG-SCAPE
conda run -n ${BigSCAPEEnv} ${bigscape} cluster \
	-i ${AntiSMASHOut} \
	-o ${BigSCAPEOut} \
	-p ${PfamDB} \
	-m 4.0 \
	--classify category \
	--alignment-mode auto \
	--extend-strategy legacy \
	--record-type region \
	--include-singletons \
	-c ${BigSCAPECpus}
