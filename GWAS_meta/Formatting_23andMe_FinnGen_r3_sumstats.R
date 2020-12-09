####################################

## Formatting summary statistics for IVW meta-analysis

## William Reay (2020)

####################################

library(dplyr)
library(data.table)
library(purrr)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/INPUT_META/")

## Common variants - default INFO (min r2 > 0.3 or mean r2 > 0.5) - 23andMe ###

##############

## Read in default info summary statistics

Default_common_23andMe <- fread("Default_INFO_23andMe_pneumonia.tab.txt.gz", header = T)

## Read in genotype info 23andMe

gt_23andMe_variants <- fread("gt_snp_stat-4.1.txt", header = T)

## Read in imputation info 23andMe

imp_23andMe_variants <- fread("im_snp_stat-4.1.txt.gz", header = T)

## Extract only variants with effect allele freq > 0.01 or < 0.99

Common_gt_var <- gt_23andMe_variants %>% filter(freq.b > 0.01 & freq.b < 0.99)
Common_imp_var <- imp_23andMe_variants %>% filter(freq.b > 0.01 & freq.b < 0.99)

Common_gt_var <- Common_gt_var %>% select(assay.name)
Common_imp_var <- Common_imp_var %>% select(assay.name)

Common_23andMe_var <- rbind(Common_gt_var, Common_imp_var)

## Merge common variants with full 23andMe summary statistics

Merged_common_default_23andMe <- merge(Common_23andMe_var, Default_common_23andMe, by="assay.name")

## Export common variants in 23andMe

write.table(Merged_common_default_23andMe, file="Common_variant_meta_input/23andMe_PNEUMONIA_common_var.tab.txt",
            sep="\t", quote = F, row.names = F)

##################

## Rare 23andMe variants - min r2 > 0.5 or mean r2 > 0.7

## Extract only variants with effect allele freq < 0.01 or > 0.99

Rare_gt_var <- gt_23andMe_variants %>% filter(freq.b < 0.01 | freq.b > 0.99)
Rare_imp_var <- imp_23andMe_variants %>% filter(freq.b < 0.01 | freq.b > 0.99)

## Filter by INFO (min > 0.5 or mean > 0.7)

Rare_imp_var <- Rare_imp_var %>% filter(min.rsqr > 0.5 | avg.rsqr > 0.7)

Rare_gt_var <- Rare_gt_var %>% select(assay.name)
Rare_imp_var <- Rare_imp_var %>% select(assay.name)

Rare_23andMe_var <- rbind(Rare_gt_var, Rare_imp_var)

## Merge rare variants with 23andMe default variants

Merged_23andMe_rare_var <- merge(Default_common_23andMe, Rare_23andMe_var, by="assay.name")

write.table(Merged_23andMe_rare_var, file="Rare_variant_meta_input/23andMe_rare_variant_input.txt", sep = "\t", row.names = F,
            quote = F)


##############

## Strict INFO for munging, min R2 > 0.9

## Filter such that min R2 > 0.9

Common_imp_var_STRICT <- Common_imp_var %>% filter(min.rsqr > 0.9)

Common_gt_var <- Common_gt_var %>% select(assay.name)
Common_imp_var_STRICT <- Common_imp_var_STRICT %>% select(assay.name)

Common_23andMe_var_STRICT <- rbind(Common_gt_var, Common_imp_var_STRICT)

Merged_common_STRICT_23andMe <- merge(Common_23andMe_var_STRICT, Default_common_23andMe, by="assay.name")

## Export common variants in 23andMe

write.table(Merged_common_STRICT_23andMe, file="INFO_over_0.9_meta_input/STRICT_INFO_23andMe_PNEUMONIA_common_var.tab.txt",
            sep="\t", quote = F, row.names = F)

################

## FINNGEN ##

## Common variants - FinnGen - mean info > 0.5 or min info > 0.3 ##

################

## Import FinnGen all pneumoniae r3 summary stats (NO INFO filtering applied yet)

FinnGen_raw <- fread("FinnGen_all_pneumoniae_r3.tab.txt", header = T)

## Import INFO file - mean INFO and per cohort INFO

FinnGen_info <- fread("FinnGen_INFO_r3.txt", header=T)

## Determine mininum INFO score for each variant

FinnGen_info <- FinnGen_info %>%
  mutate(Min = select(., starts_with("INFO")) %>% reduce(pmin))

## Retain variants with a min INFO > 0.3 OR a mean INFO > 0.5

INFO_filtered_FinnGen <- FinnGen_info %>% filter(Min > 0.3 | INFO > 0.5)

## Rename var column to match raw summary stats and merge

INFO_filtered_FinnGen <- rename(INFO_filtered_FinnGen, "chrom:pos:ref:alt"="variant")

INFO_IDs <- INFO_filtered_FinnGen %>% select(`chrom:pos:ref:alt`)

Merged_INFO_filtered_summstats <- merge(INFO_IDs, FinnGen_raw, by="chrom:pos:ref:alt")

## Filter MAF < 0.01 or MAF > 0.99

Common_FinnGen <- Merged_INFO_filtered_summstats %>% filter(maf > 0.01 & maf < 0.99)

write.table(Common_FinnGen, file="Common_variant_meta_input/FinnGen_common_var_info_filtered.tab.txt",
            sep = "\t", quote = F, row.names = F)

#############################

## Rare variants min info > 0.5 or mean info > 0.7

Rare_INFO_filtered_FinnGen <- FinnGen_info %>% filter(Min > 0.5 | INFO > 0.7)


Rare_INFO_filtered_FinnGen <- rename(Rare_INFO_filtered_FinnGen, "chrom:pos:ref:alt"="variant")

Rare_INFO_IDs <- Rare_INFO_filtered_FinnGen %>% select(`chrom:pos:ref:alt`)

Merged_rare_INFO_filtered_summstats <- merge(Rare_INFO_IDs, FinnGen_raw, by="chrom:pos:ref:alt")

## Filter MAF > 0.01 or MAF < 0.99

Rare_FinnGen <- Merged_rare_INFO_filtered_summstats %>% filter(maf < 0.01 | maf > 0.99)

write.table(Rare_FinnGen, file="Rare_variant_meta_input/FinnGen_rare_variants.tab.txt",
            sep = "\t", quote = F, row.names = F)

#######################

## FinnGen - strict for munging - min r2 > 0.9 ##

STRICT_INFO_filtered_FinnGen <- FinnGen_info %>% filter(Min > 0.9)

STRICT_INFO_filtered_FinnGen <- rename(STRICT_INFO_filtered_FinnGen, "chrom:pos:ref:alt"="variant")

STRICT_INFO_IDs <- STRICT_INFO_filtered_FinnGen %>% select(`chrom:pos:ref:alt`)

Merged_STRICT_INFO_filtered_summstats <- merge(STRICT_INFO_IDs, FinnGen_raw, by="chrom:pos:ref:alt")

## Filter MAF < 0.01 or MAF > 0.99

Common_FinnGen_STRICT <- Merged_STRICT_INFO_filtered_summstats %>% filter(maf > 0.01 & maf < 0.99)

write.table(Common_FinnGen_STRICT, file="INFO_over_0.9_meta_input/FinnGen_STRICT_INFO.tab.txt",
            sep = "\t", quote = F, row.names = F)
