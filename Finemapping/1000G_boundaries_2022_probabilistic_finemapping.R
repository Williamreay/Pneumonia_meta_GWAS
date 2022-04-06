##############################

## 23andMe/FinnGen meta 2022 

## Probabilistic finemapping - 1000 genomes LD ref locus boundaries

## William Reay (2022)

##############################

library(dplyr)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/")

## Read in ABF function

source("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/Scripts/Pneumonia_meta_GWAS_drug_repurposing/Finemapping/ABF_function.R")

## Read in sumstats

Raw_sumstats <- fread("Common_var_results/CHR_BP_annotated_FINAL_IVW_meta.txt.gz")

## Loci - chr1:154395212-154428283
##      - chr5:40486896-40524860
##      - chr11:1110395-1225078
##      - chr12:6440009-6455098

## Chromosome 1

Chr_1_subset <- Raw_sumstats %>% filter(CHR == "chr1" & BP > 154395212 & BP < 154428283)

Chr1_finemap <- Finemapping_abf(Chr_1_subset$rsID, Chr_1_subset$Effect, Chr_1_subset$StdErr)

write.table(Chr1_finemap, file="Common_var_results/Finemapping/1000_genomes_based_coord/1000G_Chr1_locus_95_percent_credible.txt",
            sep = "\t", row.names = F, quote = F)

## Chromosome 5

Chr_5_subset <- Raw_sumstats %>% filter(CHR == "chr5" & BP > 40486896 & BP < 40524860)

Chr5_finemap <- Finemapping_abf(Chr_5_subset$rsID, Chr_5_subset$Effect, Chr_5_subset$StdErr)

write.table(Chr5_finemap, file="Common_var_results/Finemapping/1000_genomes_based_coord/1000G_Chr5_locus_95_percent_credible.txt",
            sep = "\t", row.names = F, quote = F)

## Chromosome 11 and 12 have same boundaries
