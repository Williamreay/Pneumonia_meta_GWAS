#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=13:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Pneumonia - GTEx v7 Lung weights

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats IVW_META_munged_TWAS_input.sumstats.gz \
--weights  ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/Lung.P01.pos \
--weights_dir ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/ \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ~/data/users/william/23andMe_pneumonia/IVW_23andMe_FinnGen_FINAL/TWAS/Lung/Pneumonia_GTEx_V7_Lung_chr${chr}_FUSION
done
