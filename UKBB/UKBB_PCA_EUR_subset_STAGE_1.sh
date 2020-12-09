#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant filtering for PCA calculation in the EUR subset

## MAF > 0.05
## Physically genotyped on both array types
## Remove regions of long range LD

for chr in {1..22};
do
  plink2 --bgen "/home/control/ukbb2/Variant_QC_FINAL/FINAL_QC_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen.bgen" ref-first \
  --extract /home/control/data/UKBB/UKBB_variant_QC/UKBB_variant_QC/SNPs_on_both_arrays.txt \
  --exclude /home/control/data/UKBB/PCA_EUR_subset/long-range-ld.txt \
  --maf 0.05 \
  --make-bed --keep-allele-order \
  --out "/home/control/ukbb2/PC_EUR_subset_calculation/PCA_input_EUR_UKBB_non_LD_pruned_chr${chr}" ;
done
