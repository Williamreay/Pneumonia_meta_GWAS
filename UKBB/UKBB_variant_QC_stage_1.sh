#!/bin/bash

#SBATCH --cpus-per-task=6
#SBATCH --time=120:00:00
#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

## Variant QC script - initial imputation and sample exclusion

## Retain HRC sites only (AF > 1e-4)
## INFO > 0.8 for common variants and INFO > 0.8 for rare variants (AF > 1e-4)
## Unrelated white british individuals who passed sample QC retained

for chr in {1..22}
do
	plink2 --bgen "/home/control/ukbb/Raw/ukb_imp_chr${chr}_v3.bgen" ref-first \
	--extract /home/control/data/UKBB/UKBB_variant_QC/HRC_INFO_filtered_sites_to_keep.tab \
	--keep White_british_unrelated_IDs_to_keep.txt \
	--sample /home/control/ukbb/ukb58432_imp_chr1_v3_s487283.sample \
	--export bgen-1.2 \
	--out "/home/control/ukbb2/Variant_QC_Stage_1/Unrelated_white_british_INFO_HRC_filtered_chr${chr}_ukb_imp.bgen" ;
done
