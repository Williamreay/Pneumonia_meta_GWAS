###################################

## mtCOJO conditional GWAS of pneumonia - Cigs per day

#####################################

library(dplyr)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/Common_var/COJO_conditional/")

## Read pneumonia data

Pneumonia <- fread("../IVW_FINAL_Pneumonia_FinnGen_23andMe.tab.txt.gz", header = T)

Pneumonia$Allele1 <- toupper(Pneumonia$Allele1)
Pneumonia$Allele2 <- toupper(Pneumonia$Allele2)

CigsperDay <- fread("CigarettesPerDay.txt.gz", header = T)

## Rename Pneumonia columns to match cigs per day

Pneumonia <- rename(Pneumonia, "RSID" = "MarkerName", "ALT" = "Allele1", "REF"="Allele2")

## Merge by SNP

Merged_Pneumonia_SmkInit <- merge(Pneumonia, CigsperDay, by="RSID")

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

Pneumonia_mtCOJO$N <- 266277

Pneumonia_mtCOJO <- rename(Pneumonia_mtCOJO, "SNP" = "RSID", "A1" = "ALT.x",
                           "A2"="REF.x", "b" = "Effect", "se" = "StdErr",
                           "p"="P-value", "freq" = "AF")

## Write output

write.table(Pneumonia_mtCOJO, file="mtCOJO_pneumonia_meta.txt", sep = "\t", row.names = F, quote = F)


SmkInit_mtCOJO <- Merged_Pneumonia_SmkInit %>% select(RSID, ALT.y, REF.y, AF_Cigs, BETA, SE, PVALUE, N)


SmkInit_mtCOJO <- rename(SmkInit_mtCOJO, "SNP" = "RSID", "A1" = "ALT.y",
                         "A2"="REF.y", "b" = "BETA", "se" = "SE",
                         "p"="PVALUE", "freq" = "AF_Cigs")

write.table(SmkInit_mtCOJO, file="mtCOJO_CigsPerDay_meta.txt", sep = "\t", row.names = F, quote = F)
