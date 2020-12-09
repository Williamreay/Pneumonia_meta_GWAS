#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=130:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant QC - STAGE 2 - FINAL

## Exclude individuals with > 2% missingness of any chromosome (N=30)

for chr in {1..22};
do
  	plink2 --bgen "/home/control/ukbb2/Variant_QC_FINAL/TEMP_FINAL/FINAL_QC_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen.bgen" ref-first \
    --sample "/home/control/ukbb2/Variant_QC_FINAL/TEMP_FINAL/FINAL_QC_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen.sample" \
    --remove /home/control/ukbb2/Variant_QC_FINAL/TEMP_FINAL/Individuals_to_remove_per_chromosome_missingness_over_2_percent.txt \
    --export bgen-1.2 \
  	--out "/home/control/ukbb2/Variant_QC_FINAL/FINAL_QC_missingness_filtered_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen" ;
done
