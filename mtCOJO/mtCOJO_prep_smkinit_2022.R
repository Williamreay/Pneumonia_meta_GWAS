###################################

## mtCOJO conditional GWAS of pneumonia - EVER SMOKED REGULARLY

#####################################

library(dplyr)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/")

## Read pneumonia data

Pneumonia <- fread("Meta_results/FINAL_Common_IVW_FinnGen_23andMe_2022.txt.gz", header = T)

Pneumonia$Allele1 <- toupper(Pneumonia$Allele1)
Pneumonia$Allele2 <- toupper(Pneumonia$Allele2)

SmkInit <- fread("../FINAL_IVW_meta_analysis/Common_var/COJO_conditional/SmokingInitiation.txt.gz", header = T)

## Rename Pneumonia columns to match cigs per day

Pneumonia <- rename(Pneumonia, "RSID" = "MarkerName", "ALT" = "Allele1", "REF"="Allele2")

## Merge by SNP

Merged_Pneumonia_SmkInit <- merge(Pneumonia, SmkInit, by="RSID")

Merged_Pneumonia_SmkInit$ALT.y <- as.character(Merged_Pneumonia_SmkInit$ALT.y)
Merged_Pneumonia_SmkInit$ALT.x <- as.character(Merged_Pneumonia_SmkInit$ALT.x)
Merged_Pneumonia_SmkInit$REF.x <- as.character(Merged_Pneumonia_SmkInit$REF.x)
Merged_Pneumonia_SmkInit$REF.y <- as.character(Merged_Pneumonia_SmkInit$REF.y)

## Save AF relative to CigsPerDay effect allele

Merged_Pneumonia_SmkInit$AF_Cigs <- Merged_Pneumonia_SmkInit$AF

## Flip pneumonia alleles to match cigarettes per day where neccessary

mismatch <- which(Merged_Pneumonia_SmkInit$ALT.y!=Merged_Pneumonia_SmkInit$ALT.x,arr.ind=TRUE)
Merged_Pneumonia_SmkInit[mismatch,]$Effect <- Merged_Pneumonia_SmkInit[mismatch,]$Effect*-1
Merged_Pneumonia_SmkInit[mismatch,]$ALT.x <- Merged_Pneumonia_SmkInit[mismatch,]$ALT.y
Merged_Pneumonia_SmkInit[mismatch,]$REF.x <- Merged_Pneumonia_SmkInit[mismatch,]$REF.y

## Format for mtCOJO [SNP, A1, A2, freq, b, se, p, N]

Pneumonia_mtCOJO <- Merged_Pneumonia_SmkInit %>% select(RSID, ALT.x, REF.x, AF, Effect, StdErr, `P-value`)

Pneumonia_mtCOJO$N <- 391044

Pneumonia_mtCOJO <- rename(Pneumonia_mtCOJO, "SNP" = "RSID", "A1" = "ALT.x",
                           "A2"="REF.x", "b" = "Effect", "se" = "StdErr",
                           "p"="P-value", "freq" = "AF")

## Write output

write.table(Pneumonia_mtCOJO, file="Common_var_results/mtCOJO/mtCOJO_pneumonia_meta.txt", sep = "\t", row.names = F, quote = F)


SmkInit_mtCOJO <- Merged_Pneumonia_SmkInit %>% select(RSID, ALT.y, REF.y, AF_Cigs, BETA, SE, PVALUE, N)


SmkInit_mtCOJO <- rename(SmkInit_mtCOJO, "SNP" = "RSID", "A1" = "ALT.y",
                         "A2"="REF.y", "b" = "BETA", "se" = "SE",
                         "p"="PVALUE", "freq" = "AF_Cigs")

write.table(SmkInit_mtCOJO, file="Common_var_results/mtCOJO/mtCOJO_SmkInit_meta.txt", sep = "\t", row.names = F, quote = F)
