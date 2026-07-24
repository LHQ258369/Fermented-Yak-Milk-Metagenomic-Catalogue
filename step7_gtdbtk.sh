#!/bin/bash

### Step7 Taxonomic classification and phylogenomic inference by GTDB-Tk

# Conda environment
GTDBTkEnv=gtdb-tk2.4.0

# Software
gtdbtk=gtdbtk

# Directories
GenomeDir=dereplicated_genomes
ClassifyOut=gtdbtk_classify
InferOut=gtdbtk_infer
GTDBTK_DATA_PATH=gtdbtk_data/release220

# Files
MSAFile=${ClassifyOut}/align/gtdbtk.bac120.msa.filt.fasta
Prefix=901-MAGs

source activate ${GTDBTkEnv}
export GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}

### Taxonomic classification by GTDB-Tk
time -p ${gtdbtk} classify_wf \
	--genome_dir ${GenomeDir} \
	--out_dir ${ClassifyOut} \
	--cpus 20 \
	--extension fa \
	--skip_ani_screen

### Phylogenomic inference by GTDB-Tk
time -p ${gtdbtk} infer \
	--msa_file ${MSAFile} \
	--out_dir ${InferOut} \
	--prefix ${Prefix} \
	--cpus 20
