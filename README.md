Title: A metagenomic atlas reveals interaction diversity of simply assembled microbiomes of fermented yak milk in the Tibetan Plateau
 
Project Overview: This repository contains a complete set of shell scripts for shotgun metagenomic analysis of traditionally fermented yak milk from the Tibet Plateau. The pipeline covers raw read quality control, host sequence decontamination, metagenomic assembly, metagenome-assembled genome (MAG) binning, dereplication, taxonomic annotation, gene quantification, secondary metabolite biosynthetic gene cluster prediction, and microbial metabolic interaction inference.  

Script List:
Step1_fastp.sh # Raw sequencing read quality trimming
Step2_remove_host.sh # Mapping and removal of host contaminant sequences
Step3_megahit.sh # Metagenomic contig assembly 
Step4_metaphlan_humann.sh # Community profiling and functional pathway quantification
Step5_binning.sh # Contig binning, MAG recovery and quality evaluation
Step6_drep.sh # Dereplication of MAGs
Step7_gtdbtk.sh # Taxonomic classification of MAGs
Step8_salmon.sh # Gene catalog quantification
Step9_antismash_bigscape.sh # Secondary metabolite BGC prediction and comparative analysis 
Step10_metabolic_interaction.sh # Step10: Prediction of microbial interactions based on metabolic networks

Execution notes: Run scripts sequentially from step1 to step10, as downstream steps depend on outputs from previous steps; Before execution, edit each script to customize working directories, database paths, and sample information.

Citation: If you have used this script in your research, please use the following link for references to our script: https://github.com/LHQ258369/Fermented-Yak-Milk-Metagenomic-Atlas.

Limitations: This workflow was designed specifically for "A metagenomic atlas reveals interaction diversity of simply assembled microbiomes of fermented yak milk in the Tibetan Plateau"; editing and revisions might be required before applying to other projects.
