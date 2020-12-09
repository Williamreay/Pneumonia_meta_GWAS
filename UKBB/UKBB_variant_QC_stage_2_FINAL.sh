#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant QC - STAGE 2 - FINAL

## Exclude variants that satisfy any of the following
## HWE P-val < 1e-10
## Missingness > 2%
## AF < 1e-4

for chr in {1..22};
do
  	plink2 --bgen "/home/control/ukbb2/Variant_QC_Stage_1/Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen" ref-first \
    --sample "/home/control/ukbb2/Variant_QC_Stage_1/Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.sample" \
    --geno 0.02 \
    --maf 0.0001 \
    --hwe 0.0000000001 \
    --export bgen-1.2 \
  	--out "/home/control/ukbb2/Variant_QC_FINAL/FINAL_QC_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen" ;
done
