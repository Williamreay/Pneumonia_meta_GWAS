#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=13:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#Pneumonia - GTEx v7 Whole_Spleen weights

for chr in $(seq 1 22);
do
Rscript ~/data/users/william/fusion_twas-master/FUSION.assoc_test.R \
--sumstats Full_munged/Munged_meta_analysis_ALL.sumstats.gz  \
--weights  ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/Spleen.P01.pos \
--weights_dir ~/data/users/william/Cytokine_pneumonia_lung_function/TWAS/TWAS_weights/GTEx.Spleen.P01 \
--ref_ld_chr ~/data/users/william/fusion_twas-master/LDREF/1000G.EUR. \
--chr ${chr} \
--out ~/data/users/william/2022_FinnGen_R6_meta/Meta_analysis_output/TWAS/2022_Spleen/2022_Pneumonia_GTEx_V7_Spleen_chr${chr}_FUSION
done
