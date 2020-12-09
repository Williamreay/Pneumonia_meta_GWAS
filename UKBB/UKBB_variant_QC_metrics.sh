#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=80:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant QC script - output variant metrics, HWE, missinginess, AF etc.

for chr in {1..22};
do
  qctool -g ~/ukbb/Variant_QC_Stage_1/Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen \
  -threads 6 \
  -snp-stats -osnp UKBB_white_british_unrelated_variant_metrics_chr${chr}.txt ;
done
