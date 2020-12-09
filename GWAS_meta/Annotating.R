library(data.table)

library(dplyr)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/")

CHR_hg19_23andMe <- fread("INPUT_META/Common_variant_meta_input/23andMe_PNEUMONIA_common_var.tab.txt", header = T)

CHR_hg19_23andMe <- CHR_hg19_23andMe %>% select(assay.name, scaffold, position, OTHER_ALLELE, EFFECT_ALLELE)

CHR_hg19_23andMe <- rename(CHR_hg19_23andMe, "MarkerName"="assay.name")

## Make input for finemapping pipeline

Common_var_IVW <- fread("Common_var/IVW_FINAL_Pneumonia_FinnGen_23andMe.tab.txt.gz", header = T)

CHR_BP_annotated <- merge(CHR_hg19_23andMe, Common_var_IVW, by="MarkerName")

CHR_BP_annotated <- unique(CHR_BP_annotated)

## Load 23andMe allele freq

## Read in genotype info 23andMe

gt_23andMe_variants <- fread("INPUT_META/gt_snp_stat-4.1.txt", header = T)

## Read in imputation info 23andMe

imp_23andMe_variants <- fread("INPUT_META/im_snp_stat-4.1.txt.gz", header = T)

gt_23andMe_variants <- gt_23andMe_variants %>% select(assay.name, freq.a, freq.b)
imp_23andMe_variants <- imp_23andMe_variants %>% select(assay.name, freq.a, freq.b)

Freq_23andMe_var <- rbind(gt_23andMe_variants, imp_23andMe_variants)

Freq_23andMe_var <- Freq_23andMe_var[!duplicated(Freq_23andMe_var[, "assay.name"])]

Freq_23andMe_var <- rename(Freq_23andMe_var, "MarkerName"="assay.name")

Freq_annotated_23andMe_raw <- merge(CHR_BP_annotated, Freq_23andMe_var, by="MarkerName")

Freq_annotated_23andMe_raw <- unique(Freq_annotated_23andMe_raw)

Freq_annotated_23andMe_raw$Allele1 <- toupper(Freq_annotated_23andMe_raw$Allele1)
Freq_annotated_23andMe_raw$Allele2 <- toupper(Freq_annotated_23andMe_raw$Allele2)

## Match Allele1 to EFFECT ALLELE (freq b) in 23andMe

mismatch_1 <- which(Freq_annotated_23andMe_raw$Allele1!=Freq_annotated_23andMe_raw$EFFECT_ALLELE,arr.ind=TRUE)
Freq_annotated_23andMe_raw[mismatch_1,]$Effect <- Freq_annotated_23andMe_raw[mismatch_1,]$Effect*-1
Freq_annotated_23andMe_raw[mismatch_1,]$Allele1 <- Freq_annotated_23andMe_raw[mismatch_1,]$EFFECT_ALLELE
Freq_annotated_23andMe_raw[mismatch_1,]$Allele2 <- Freq_annotated_23andMe_raw[mismatch_1,]$OTHER_ALLELE

## Flip effect allele if freq.b > 0.5

mismatch_2 <- which(Freq_annotated_23andMe_raw$freq.b > 0.5, arr.ind = TRUE)
Freq_annotated_23andMe_raw[mismatch_2,]$Effect <- Freq_annotated_23andMe_raw[mismatch_2,]$Effect*-1
Freq_annotated_23andMe_raw[mismatch_2,]$Allele1 <- Freq_annotated_23andMe_raw[mismatch_2,]$OTHER_ALLELE
Freq_annotated_23andMe_raw[mismatch_2,]$Allele2 <- Freq_annotated_23andMe_raw[mismatch_2,]$EFFECT_ALLELE
Freq_annotated_23andMe_raw[mismatch_2,]$freq.b <- Freq_annotated_23andMe_raw[mismatch_2,]$freq.a

## CHR, BP, rsID, MAF (0-0.5), EA, NEA, BETA, SE, P, 

Freq_annotated_23andMe_raw <- rename(Freq_annotated_23andMe_raw, "CHR"="scaffold", "BP"="position", "rsID"="MarkerName",
                                   "EA"="Allele1", "NEA"="Allele2", "P"="P-value",  "MAF"="freq.b")

                               
Freq_annotated_23andMe_raw <- Freq_annotated_23andMe_raw %>% select(CHR, BP, rsID, MAF, EA, NEA, Effect, StdErr, P)


## write output for all SNPs

write.table(Freq_annotated_23andMe_raw, file="~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/Common_var/echolocatoR_pipeline_input/Finemapping_IVW_common_var_pneumonia_input.txt",
            sep="\t", row.names = F, quote = F)
