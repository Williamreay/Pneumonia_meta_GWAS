################################

## Preparing phenotypes for UKBB PRS/PES analyses

## William Reay

#################################

library(dplyr)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Pneumonia_phenotypes/")

source("Pneumonia_phenotypes.r")

## Load data 

## Self-reported or ICD-10

Self_reported_or_ICD10 <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Self_reported_OR_ICD10_pneumonia.tab.txt", header = T)

## ICD-10 (self-reported removed as controls)

ICD10 <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/ICD10_pneumonia_self_reported_pneumonia_removed_as_controls.tab.txt", header = T)

## IDs in genetics

Gen_IDs <- fread("~/Desktop/UKBB_info/Ethinicity_relatedness_QC/White_british_unrelated_IDs_to_keep.txt", header = T)


## Merge

Merged_self <- merge(Self_reported_or_ICD10, Gen_IDs, by="f.eid")

Merged_self <- rename(Merged_self, "IID"="f.eid")

write.table(Merged_self, file="~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Self_reported_or_ICD10_in_genetics.tab.txt",
            sep="\t", row.names = F, quote = F)

Merged_ICD <- merge(ICD10, Gen_IDs, by="f.eid")

Merged_ICD <- rename(Merged_ICD, "IID"="f.eid")

write.table(Merged_ICD, file="~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/ICD10_in_genetics.tab.txt",
            sep="\t", row.names = F, quote = F)





