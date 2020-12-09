#!/bin/bash

#SBATCH --cpus-per-task=6
#SBATCH --time=120:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## PCA with flashPCA2 on the filtered genotypes as per the previous three steps

flashpca --bfile ~/ukbb2/PC_EUR_subset_calculation/LD_pruned/Merged_PCA_input_EUR_UKBB_LD_pruned -d 20 -v 
