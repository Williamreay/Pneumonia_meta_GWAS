########################################

## CRP ==> pneumonia susceptibility (2022)

## Non-UKBB CRP GWAS

## Univariable MR - total effect

## William Reay (2022)

########################################


library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/LCV/MR_CRP_follow_up/")

## Get non-UKBB CRP IVs

CRP_IVs <- extract_instruments("ieu-b-35", clump = TRUE)

Pneumonia_raw <- fread("../../../Meta_results/FINAL_Common_IVW_FinnGen_23andMe_2022.txt.gz")
Pneumonia_raw$Phenotype <- "Pneumonia_susceptibility"

Pneumonia_outcome <-  format_data(dat = Pneumonia_raw,
                                  type = "outcome",
                                  beta_col = "Effect",
                                  se_col = "StdErr",
                                  snp_col = "MarkerName",
                                  effect_allele_col = "Allele1",
                                  other_allele_col = "Allele2",
                                  phenotype_col = "Phenotype")


## Harmonise data and exclude potential palindromes (option 3)

CRP_pneu_harm <- harmonise_data(CRP_IVs, Pneumonia_outcome, action = 3)

Univariable_CRP <- mr(CRP_pneu_harm,
                      method_list = c("mr_ivw_mre",
                                      "mr_ivw_fe",
                                      "mr_weighted_median",
                                      "mr_weighted_mode",
                                      "mr_egger_regression"))

OR_CRP_to_pneumonia <- generate_odds_ratios(Univariable_CRP)

write.table(OR_CRP_to_pneumonia, file="Non_UKBB_Univariable_CRP_to_pneumonia.txt",
            sep = "\t", row.names = F, quote = F)

## Heterogeneity test via Cochran's Q

CRP_pneumonia_het <- mr_heterogeneity(CRP_pneu_harm)

## Test if MR-Egger intercept is significantly different from zero

CRP_pneumonia_Egger_intercept <- mr_pleiotropy_test(CRP_pneu_harm)

CRP_pneumonia_LOO_IVW <- mr_leaveoneout(CRP_pneu_harm, method = mr_ivw_mre)

CRP_pneumonia_LOO_IVW %>% filter(p > 0.05)

CRP_single_SNP <- mr_singlesnp(CRP_pneu_harm)

## Derive top three single SNPs

CRP_single_SNP <-  CRP_single_SNP[order(CRP_single_SNP$p),]

mr_forest_plot(CRP_single_SNP, exponentiate = F)
