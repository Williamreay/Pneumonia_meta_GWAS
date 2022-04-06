#!/bin/bash

#PBS -l select=1:ncpus=16:mem=150G,walltime=200:00:00

## PRSet - MUC5AC PES permutation test

/home/c3202333/PRSice_linux --A1 EA --A2 NEA \
  --chr CHR --bp BP --stat Effect \
  --snp rsID --thread 15 \
  --base /home/c3202333/Pneumonia_PRS_PES/MUC5AC_weights/Scoring_input_pneumonia_GWAS.txt \
  --target /home/c3202333/Variant_QC_FINAL/Biallelic_common_merged_plink2_binary/COMMON_BIALLELIC_plink_bed \
  --binary-target T \
  --pheno /home/c3202333/Pneumonia_PRS_PES/MUC5AC_weights/Pneumonia_BROAD_phenotype.pheno \
  --beta \
  --fastscore \
  --all-score \
  --bar-levels 0.005 \
  --no-full \
  --gtf /home/c3202333/Pneumonia_PRS_PES/Homo_sapiens.GRCh37.87.gtf \
  --msigdb /home/c3202333/Pneumonia_PRS_PES/MsigDB_canonical_and_hallmark_Tclin_v_6.1.0.gmt \
  --set-perm 10000 --wind-3 35000 --wind-5 10000 \
  --out /home/c3202333/Pneumonia_PRS_PES/Tclin_all_MSigDB/Liberal_boundaries_PRSET_0_005 \
  --no-full 
