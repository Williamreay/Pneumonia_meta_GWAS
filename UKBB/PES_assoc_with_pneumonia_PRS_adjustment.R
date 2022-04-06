#############################################

## Pneumonia PES - MUC5AC overlapping

## William Reay (2022)

# Broad phenotype

#############################################


library(dplyr)
library(data.table)
library(rcompanion)
library(caret)
library(easyGgplot2)


Covariates <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)

## Load broad phenotype definition

Broad_pneumonia <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Self_reported_or_ICD10_in_genetics.tab.txt")

Merged_Broad <- merge(Covariates, Broad_pneumonia, by = "IID")

Merged_Broad$Batch <- as.factor(Merged_Broad$Batch)

## PES

PES_raw <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/Tclin_UKBB_PRset_permutation/Scores_for_immune_corr/MUC5AC_overlapping_replicate.txt")

Colnames_PES_raw <- make.unique(names(PES_raw))

colnames(PES_raw) <- Colnames_PES_raw


## Broad

BROAD_merge <- merge(PES_raw, Merged_Broad, by = "IID")

BROAD_merge[,c(2:4)] <- lapply(BROAD_merge[,c(2:4)], function(x) c(scale(x)))

## Test adjusted for GW PRS

Protein_metab <- glm(Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                      PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch + PRS_0_05 + REACTOME_METABOLISM_OF_PROTEINS_0_05, family = "binomial",
                    data = BROAD_merge)

anova(Protein_metab, test="Chisq")

## P = 0.01

PTM <- glm(Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
             PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch + PRS_0_05 + REACTOME_POST_TRANSLATIONAL_PROTEIN_MODIFICATION_0_05, family = "binomial",
           data = BROAD_merge)

anova(PTM, test="Chisq")

## P = 0.043
