####################################

## PES - Bile acid metabolism

## William Reay (2020)

##################################

library(data.table)
library(dplyr)
library(rcompanion)
library(jtools)
library(sjPlot)
library(sjmisc)
theme_set(theme_sjplot())

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/")

## Import the bile acid metabolism PES at each P threshold on all 22 autosomes, concatenated into the sample file

All_chr_bile_PES <- fread("PES_scores/Bile_acid_metabolism/All_chr_bile_acid_metabolism.txt", sep = " ", header = T)

Colnames_all_chr_bile_PES <- make.unique(names(All_chr_bile_PES))

colnames(All_chr_bile_PES) <- Colnames_all_chr_bile_PES

## Convert double IID to single IID type

All_chr_bile_PES$IID <- gsub(".*_","",All_chr_bile_PES$IID)  

## Select PES on each chromosome

Bile_acid_PES <- All_chr_bile_PES %>% select(IID, starts_with("Pt_0.005"))

## Sum chromosome-wise PES

Bile_acid_PES$Bile_acid_PES <- rowSums(Bile_acid_PES[,-1])

Bile_acid_PES <- Bile_acid_PES %>% select(IID, Bile_acid_PES)

write.table(Bile_acid_PES, file="PES_scores/Bile_acid_metabolism/Bile_acid_metabolism_pneumonia_PES.tab.txt",
            sep = "\t", row.names = F, quote = F)

## Convert IID to numeric

Bile_acid_PES$IID <- as.numeric(Bile_acid_PES$IID)

## Test association with and without covariation for PRS

## Load covariate data

Covariates <- fread("Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)

## Load strict phenotype definition (ICD-10 only, self-reported removed as controls)

Strict_pneumonia <- fread("ICD10_in_genetics.tab.txt", header = T)

## Load broad phenotype definition (ICD-10 or self-reported as cases)

Broad_pneumonia <- fread("Self_reported_or_ICD10_in_genetics.tab.txt", header = T)

## Load raw, unscaled, PRS

PRS_all <- fread("PRS_scores/All_RAW_PRS_scores_all_chromosome.tab.txt")

## Scale PES and PRS at same threshold (P < 0.005) to have zero mean and unit variance

PRS_all$Scaled_PRS_0.005 <- as.numeric(scale(PRS_all$PRS_0.005))
Bile_acid_PES$Scaled_bile_acid_PES <- as.numeric(scale(Bile_acid_PES$Bile_acid_PES))

########### STRICT PHENOTYPE #############

## Merge covariates, pneumonia, PRS, and PES

Strict_pneumonia_phenotype_cov_PES <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                                             list(PRS_all, Covariates, Strict_pneumonia, Bile_acid_PES))

Strict_pneumonia_phenotype_cov_PES$ICD10_pneumonia <- as.factor(Strict_pneumonia_phenotype_cov_PES$ICD10_pneumonia)

## Null model (intercept + covariates)

Strict_NULL <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     Batch,family="binomial", data = Strict_pneumonia_phenotype_cov_PES)

## Test association with PES by itself (covaried for same variables as PRS models)

Strict_Bile_no_PRS <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                     Batch + Scaled_bile_acid_PES,family="binomial", data = Strict_pneumonia_phenotype_cov_PES)

## Test association of top decile vs rest of cohort

Strict_pneumonia_phenotype_cov_PES$Decile_PES <- ntile(Strict_pneumonia_phenotype_cov_PES$Bile_acid_PES,10)

Strict_pneumonia_phenotype_cov_PES$Top_decile_PES <- ifelse(Strict_pneumonia_phenotype_cov_PES$Decile_PES == 10, 1, 0)

Decile_strict_Bile_no_PRS <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            Batch + Top_decile_PES,family="binomial", data = Strict_pneumonia_phenotype_cov_PES)


########### BROAD PHENOTYPE #############

## Merge covariates, pneumonia, and PRS

Broad_pneumonia_phenotype_cov_PES <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                                            list(PRS_all, Covariates, Broad_pneumonia, Bile_acid_PES))

Broad_pneumonia_phenotype_cov_PES$Self_reported_or_ICD <- as.factor(Broad_pneumonia_phenotype_cov_PES$Self_reported_or_ICD)

## Create null model with intercept + covariates only

Broad_NULL <- glm(Self_reported_or_ICD ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Batch,
                  family="binomial", data = Broad_pneumonia_phenotype_cov_PES)

## Test association

Broad_Bile_no_PRS <- glm(Self_reported_or_ICD ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 +
                            Batch + Scaled_bile_acid_PES,family="binomial", data = Broad_pneumonia_phenotype_cov_PES)


## Correlation with PES and PRS at the same P threshold

PES_PRS_lm <- cor(Scaled_PRS_0.005 ~ Scaled_bile_acid_PES, data = Broad_pneumonia_phenotype_cov_PES)

## Test whether for individuals in the top decile of bile acid metabolism PES, taking statins reduces their odds of pneumonia

## Read in individuals who self-report taking statins

Statins <- fread("Statins_UKBB.tab.txt", header = T)

Statins <- dplyr::rename(Statins, "IID"="f.eid")


Statins_pneumonia <- merge(Statins, Broad_pneumonia_phenotype_cov_PES,by = "IID")

## Test association between statin usage and pneumonia bile acid PES

Statin_pneumonia_test <- glm(Self_reported_or_ICD ~ Sex + Age + Age2 + Statins*Scaled_bile_acid_PES + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = Statins_pneumonia)

## Test association between PES and statin usage

Statin_PES <- glm(Statins ~ Sex + Age + Age2 + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Scaled_bile_acid_PES, family = "binomial", data = Statins_pneumonia)

## Adjust above model for PRS at same p value threshold

Statin_PES_PRS <- glm(Statins ~ Sex + Age + Age2 + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Scaled_bile_acid_PES + Scaled_PRS_0.005, family = "binomial", data = Statins_pneumonia)

exp(cbind(coef(Statin_PES_PRS), confint(Statin_PES_PRS)))

## Plot interaction between statins and PES

plot_model(Statin_pneumonia_test, type = "eff", terms = c("Scaled_bile_acid_PES", "Statins"))


## Repeat for strict phenotype


Strict_statins_pneumonia <- merge(Statins, Strict_pneumonia_phenotype_cov_PES,by = "IID")

## Test association between statin usage and pneumonia bile acid PES

Strict_statin_pneumonia_test <- glm(ICD10_pneumonia ~ Sex + Age + Age2 + Statins*Scaled_bile_acid_PES + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = Strict_statins_pneumonia)

Strict_statin_pneumonia_test_three_way <- glm(ICD10_pneumonia ~ Age + Age2 + Statins*Scaled_bile_acid_PES*Sex + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family="binomial", data = Strict_statins_pneumonia)

## Test association between PES and statin usage

Strict_satin_PES <- glm(Statins ~ Sex + Age + Age2 + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Scaled_bile_acid_PES, family = "binomial", data = Strict_statins_pneumonia)

## Adjust above model for PRS at same p value threshold

Strict_statin_PES_PRS <- glm(Statins ~ Sex + Age + Age2 + Batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Scaled_bile_acid_PES + Scaled_PRS_0.005, family = "binomial", data = Strict_statins_pneumonia)


plot_model(Strict_statin_pneumonia_test, type = "eff", terms = c("Scaled_bile_acid_PES", "Statins"))

plot_model(Strict_statin_pneumonia_test_three_way, type = "eff", terms = c("Scaled_bile_acid_PES", "Statins", "Sex"))
