#!/bin/bash

### Step6 Genome dereplication by dRep

# Conda environment
DrepEnv=drep

# Software
dRep=dRep

# Directories
InputBins=input
Output95=drep_95
Output99=drep_99

mkdir -p ${Output95} ${Output99}

### Activate dRep environment
source activate ${DrepEnv}

### Dereplicate bins at 95% ANI
time -p ${dRep} dereplicate ${Output95} \
	-g ${InputBins}/*.fa \
	-p 20 \
	-sa 0.95

### Dereplicate bins at 99% ANI
time -p ${dRep} dereplicate ${Output99} \
	-g ${InputBins}/*.fa \
	-p 20 \
	-sa 0.99
