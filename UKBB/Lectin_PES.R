####################################

## PES - Lectin pathway

## William Reay (2020)

##################################

library(data.table)
library(dplyr)
library(rcompanion)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/")

## Import the Lectin PES at each P threshold on all 22 autosomes, concatenated into the sample file

All_chr_Lectin_PES <- fread("PES_scores/Lectin_pathway/All_chr_lectin_pathway.txt", sep = " ", header = T)

Colnames_all_chr_Lectin_PES <- make.unique(names(All_chr_Lectin_PES))

colnames(All_chr_Lectin_PES) <- Colnames_all_chr_Lectin_PES

## Convert double IID to single IID type

All_chr_Lectin_PES$IID <- gsub(".*_","",All_chr_Lectin_PES$IID)  

## Select PES on each chromosome

Lectin_PES <- All_chr_Lectin_PES %>% select(IID, starts_with("Pt_0.05"))

## Sum chromosome-wise PES

Lectin_PES$Lectin_pathway_PES <- rowSums(Lectin_PES[,-1])

Lectin_PES <- Lectin_PES %>% select(IID, Lectin_PES)

write.table(Lectin_PES, file="PES_scores/Lectin_pathway/Lectin_pathway_pneumonia_PES.tab.txt",
            sep = "\t", row.names = F, quote = F)

## Convert IID to numeric

Lectin_PES$IID <- as.numeric(Lectin_PES$IID)

## Test association with and without covariation for PRS

## Load covariate data

Covariates <- fread("Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)

## Load strict phenotype definition (ICD-10 only, self-reported removed as controls)

Strict_pneumonia <- fread("ICD10_in_genetics.tab.txt", header = T)

## Load broad phenotype definition (ICD-10 or self-reported as cases)

Broad_pneumonia <- fread("Self_reported_or_ICD10_in_genetics.tab.txt", header = T)

## Load raw, unscaled, PRS

PRS_all <- fread("PRS_scores/All_RAW_PRS_scores_all_chromosome.tab.txt")

## Scale PES and PRS at same threshold (P < 0.05) to have zero mean and unit variance

PRS_all$Scaled_PRS_0.05 <- as.numeric(scale(PRS_all$PRS_0.05))
Lectin_PES$Lectin_PES <- as.numeric(scale(Lectin_PES$Lectin_pathway_PES))

########### STRICT PHENOTYPE #############

## Merge covariates, pneumonia, PRS, and PES

Strict_pneumonia_phenotype_cov_PES <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                                             list(PRS_all, Covariates, Strict_pneumonia, Lectin_PES))

Strict_pneumonia_phenotype_cov_PES$ICD10_pneumonia <- as.factor(Strict_pneumonia_phenotype_cov_PES$ICD10_pneumonia)

## Null model (intercept + covariates)

Strict_NULL <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     Batch,family="binomial", data = Strict_pneumonia_phenotype_cov_PES)

## Test association with PES by itself (covaried for same variables as PRS models)

Strict_Lectin_no_PRS <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                  Batch + Lectin_PES,family="binomial", data = Strict_pneumonia_phenotype_cov_PES)


########### BROAD PHENOTYPE #############

## Merge covariates, pneumonia, and PRS

Broad_pneumonia_phenotype_cov_PES <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                                            list(PRS_all, Covariates, Broad_pneumonia, Lectin_PES))

Broad_pneumonia_phenotype_cov_PES$Self_reported_or_ICD <- as.factor(Broad_pneumonia_phenotype_cov_PES$Self_reported_or_ICD)


## Test association

Broad_Lectin_no_PRS <- glm(Self_reported_or_ICD ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                                 Batch + Lectin_PES,family="binomial", data = Broad_pneumonia_phenotype_cov_PES)
