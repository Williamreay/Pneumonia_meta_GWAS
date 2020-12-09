##################################
## LDSC with UKBB metabolites and pneumonia
## William Reay 2020
#################################

#!/bin/bash

for i in $(cat LDSC_input.txt);
do
~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/ldsc.py  \
--ref-ld-chr ~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/eur_w_ld_chr/ \
--w-ld-chr ~/data/users/william/Cytokine_pneumonia_lung_function/ldsc/eur_w_ld_chr/ \
--rg IVW_META_munged_HapMap3.sumstats.gz,../../SZ_metabolites_2020/IRNT_metabolites/${i} \
--out Biochem_LDSC/Pneumonia_${i}_rg_LDSC;
done
