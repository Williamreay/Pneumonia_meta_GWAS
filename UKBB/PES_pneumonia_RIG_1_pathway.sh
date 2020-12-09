#!/bin/bash

#SBATCH --cpus-per-task=6
#SBATCH --time=180:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## PRSice2 - PES: RIG-1 receptor pathway

for chr in {1..22};
do
  PRSice_linux --A1 EA --A2 NEA \
  --chr CHR --bp BP --stat Effect \
  --snp rsID \
  --base /home/control/data/users/william/23andMe_pneumonia/IVW_23andMe_FinnGen_FINAL/PES_PRS/PES_pathway_variants/RIG_PATHWAY_variant_input \
  --type bgen \
  --target "/home/control/ukbb2/Variant_QC_FINAL/FINAL_QC_Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen,/home/control/ukbb2/PRS/PRS_pneumonia_sample_file.bgen.sample" \
  --fastscore \
  --no-regress \
  --all-score \
  --no-full \
  --bar-levels 0.005 \
  --allow-inter \
  --score sum \
  --beta \
  --binary-target T \
  --model add \
  --ignore-fid \
  --print-snp \
  --extract /home/control/ukbb2/PRS/Valid_SNPs/PRS_pneumonia_self_reported_or_ICD10_chr${chr}.valid \
  --thread 11 \
  --out RIG_1_pathway/PES_RIG_1_pathway_chr${chr} ;
done
