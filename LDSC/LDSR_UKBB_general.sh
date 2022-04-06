#!/bin/bash

#SBATCH --cpus-per-task=3
#SBATCH --time=12:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Ben Neale UKBB traits LDSR

eval "$(conda shell.bash hook)"
conda activate ldsc

for i in $(cat ~/data/users/william/2022_FinnGen_R6_meta/Meta_analysis_output/LDSR_UKBB_Neale_general/LDSR_Neale_UKBB/LDSR_input_general.txt);
do
~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/ldsc.py  \
--ref-ld-chr ~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/eur_w_ld_chr/ \
--w-ld-chr ~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/eur_w_ld_chr/ \
--rg ~/data/users/william/2022_FinnGen_R6_meta/Meta_analysis_output/HapMap3_munged/Munged_meta_analysis_HapMap3.sumstats.gz,~/data/users/william/2022_FinnGen_R6_meta/Meta_analysis_output/LDSR_UKBB_Neale_general/LDSR_Neale_UKBB/${i} \
--out Pneumonia_2022_${i}_rg_UKBB_general;
done

