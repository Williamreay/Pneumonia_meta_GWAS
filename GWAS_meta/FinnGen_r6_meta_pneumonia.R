####################################

## Formatting summary statistics for IVW meta-analysis

## William Reay (2022)

####################################

library(dplyr)
library(data.table)
library(purrr)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/INPUT_META/")

## 23andMe

## Common Variants 

## INFO > 0.6

## Read in default info summary statistics

Default_common_23andMe <- fread("Default_INFO_23andMe_pneumonia.tab.txt.gz", header = T)

## Read in genotype info 23andMe

gt_23andMe_variants <- fread("gt_snp_stat-4.1.txt", header = T)

## Read in imputation info 23andMe

imp_23andMe_variants <- fread("im_snp_stat-4.1.txt.gz", header = T)

## Filter INFO > 0.6

imp_23andMe_variants <- imp_23andMe_variants %>% filter(avg.rsqr > 0.6)

## Extract only variants with effect allele freq > 0.01 or < 0.99

Common_gt_var <- gt_23andMe_variants %>% filter(freq.b > 0.01 & freq.b < 0.99)
Common_imp_var <- imp_23andMe_variants %>% filter(freq.b > 0.01 & freq.b < 0.99)

Common_gt_var <- Common_gt_var %>% select(assay.name)
Common_imp_var <- Common_imp_var %>% select(assay.name)

Common_23andMe_var <- rbind(Common_gt_var, Common_imp_var)

## Merge common variants with full 23andMe summary statistics

Merged_common_default_23andMe <- merge(Common_23andMe_var, Default_common_23andMe, by="assay.name")

Merged_common_default_23andMe <- unique(Merged_common_default_23andMe)

write.table(Merged_common_default_23andMe, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/Common_var/Common_var_23andMe_INFO_0.6.txt",
            sep = "\t", row.names = F, quote = F)

## Rare variants - INFO > 0.6

Rare_gt_var <- gt_23andMe_variants %>% filter(freq.b < 0.01 | freq.b > 0.99)
Rare_imp_var <- imp_23andMe_variants %>% filter(freq.b < 0.01 | freq.b > 0.99)
Rare_imp_var <- Rare_imp_var %>% filter(avg.rsqr > 0.6)

Rare_gt_var <- Rare_gt_var %>% select(assay.name)
Rare_imp_var <- Rare_imp_var %>% select(assay.name)

Rare_23andMe_var <- rbind(Rare_gt_var, Rare_imp_var)

## Merge

Merged_23andMe_rare_var <- merge(Default_common_23andMe, Rare_23andMe_var, by="assay.name")

write.table(Merged_23andMe_rare_var, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/Rare_var/Rare_var_23andMe_INFO_0.6.txt", sep = "\t", row.names = F,
            quote = F)

## Strict info for munging - INFO > 0.9

Strict_imp_23andMe_variants <- imp_23andMe_variants %>% filter(min.rsqr > 0.9)

Common_imp_var_STRICT <- Strict_imp_23andMe_variants %>% select(assay.name)

Common_23andMe_var_STRICT <- rbind(Common_gt_var, Common_imp_var_STRICT)

Merged_common_STRICT_23andMe <- merge(Common_23andMe_var_STRICT, Default_common_23andMe, by="assay.name")

## Export common variants in 23andMe

write.table(Merged_common_STRICT_23andMe, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/INFO_over_0.9_common_var/STRICT_23andMe_common_var_INFO_over_0.9",
            sep="\t", quote = F, row.names = F)


## FinnGen ##

## Common var - default INFO > 0.6 (already filtered)

## Import FinnGen all pneumoniae r3 summary stats (NO INFO filtering applied yet)

FinnGen_raw <- fread("../../2022_FinnGen_r6_meta_analysis/Meta_input/FinnGen_raw_pneumonia_J10.txt.gz", header = T)

## Import INFO file - mean INFO

FinnGen_info <- fread("../../2022_FinnGen_r6_meta_analysis/Meta_input/R6_FinnGen_INFO.txt", header=T)

FinnGen_info <- rename(FinnGen_info, "#chrom:pos:ref:alt"="#variant")

INFO_filtered_FinnGen <- FinnGen_info %>% filter(INFO > 0.6)


INFO_IDs <- INFO_filtered_FinnGen %>% select(`#chrom:pos:ref:alt`)

Merged_INFO_filtered_summstats <- merge(INFO_IDs, FinnGen_raw, by="#chrom:pos:ref:alt")

## Filter MAF < 0.01 or MAF > 0.99

Common_FinnGen <- Merged_INFO_filtered_summstats %>% filter(af_alt > 0.01 & af_alt < 0.99)

Common_FinnGen <- unique(Common_FinnGen)

write.table(Common_FinnGen, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/Common_var/Common_var_FinnGen_INFO_0.6.txt",
            sep = "\t", quote = F, row.names = F)

## Rare variants - default INFO > 0.6

Rare_FinnGen <- Merged_INFO_filtered_summstats %>% filter(af_alt < 0.01 | af_alt > 0.99)

Rare_FinnGen <- unique(Rare_FinnGen)

write.table(Rare_FinnGen, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/Rare_var/Rare_var_FinnGen_INFO_0.6.txt",
            sep = "\t", quote = F, row.names = F)

## Strict INFO > 0.9

Strict_INFO_filtered_FinnGen <- FinnGen_info %>% filter(INFO > 0.9)


Strict_INFO_IDs <- Strict_INFO_filtered_FinnGen %>% select(`#chrom:pos:ref:alt`)

Strict_Merged_INFO_filtered_summstats <- merge(Strict_INFO_IDs, FinnGen_raw, by="#chrom:pos:ref:alt")

## Filter MAF < 0.01 or MAF > 0.99

Strict_Common_FinnGen <- Strict_Merged_INFO_filtered_summstats %>% filter(af_alt > 0.01 & af_alt < 0.99)

Strict_Common_FinnGen <- unique(Strict_Common_FinnGen)

write.table(Strict_Common_FinnGen, file="../../2022_FinnGen_r6_meta_analysis/Meta_input/INFO_over_0.9_common_var/STRICT_FinnGen_common_var_INFO_over_0.9.txt",
            sep = "\t", quote = F, row.names = F)
