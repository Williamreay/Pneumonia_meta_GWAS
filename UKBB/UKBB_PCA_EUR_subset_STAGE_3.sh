#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant filtering for PCA calculation in the EUR subset

for chr in {1..22};
do
  plink2 --bfile "/home/control/ukbb2/PC_EUR_subset_calculation/PCA_input_EUR_UKBB_non_LD_pruned_chr${chr}" \
  --extract "/home/control/ukbb2/PC_EUR_subset_calculation/LD_pruned_chr${chr}.prune.in" \
  --make-bed \
  --out "/home/control/data/UKBB/PCA_EUR_subset/LD_pruned/PCA_input_LD_pruned_chr${chr}" ;
done
