#!/bin/bash

## GWAS of pneumonia conditioned on 'Ever smoked phenotype'

./gcta64 --bfile ~/cloudstor/Cross_disorder_PES_2019/GEUVADIS_gene_expression/g1000_eur \
--mtcojo-file smoking_init_mtcojo_summary.list \
--ref-ld-chr eur_w_ld_chr/ --w-ld-chr eur_w_ld_chr/ \
--out Pneumonia_meta_SmkInit
