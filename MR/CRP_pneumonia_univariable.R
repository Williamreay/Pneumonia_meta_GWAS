########################################

## CRP ==> pneumonia susceptibility (2022)

## Univariable MR - total effect

## William Reay (2022)

########################################

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/LCV/MR_CRP_follow_up/")

## Read in CRP exposure data in SD units

CRP_raw <- fread("CAUSE_formatted_CRP_30710_irnt.gwas.imputed_v3.both_sexes.tsv.gz")

## Retain GW sig var

CRP_gw_sig <- CRP_raw %>% filter(pval < 5e-08)

CRP_gw_sig$Phenotype <- "CRP"

## Format as IEUGWAS IV

CRP_exp <- format_data(dat = CRP_gw_sig,
                      type = "exposure",
                      beta_col = "beta",
                      se_col = "se",
                      snp_col = "rsid",
                      effect_allele_col = "alt",
                      other_allele_col = "ref",
                      phenotype_col = "Phenotype")

## LD clump per usual r2

CRP_clumped <- clump_data(CRP_exp, pop = "EUR")

## Read in pneumonia outcome data

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

CRP_pneu_harm <- harmonise_data(CRP_clumped, Pneumonia_outcome, action = 3)


## Estimate variance explained and calculate i) F-stat, and ii) I^2 stat

R2_input <- CRP_pneu_harm %>% filter(mr_keep == TRUE)

Approximated_r2_gamma_glutamyltransferase <- sum(get_r_from_pn(R2_input$pval.exposure,343524)^2)

## Calculate F statistic to evaluate whether it's > 10

F_statistic <- function(r,n,k){r*(n-1-k)/(1-r)/k}

F_statistic(Approximated_r2_gamma_glutamyltransferase, 343524, 156)

## F stat = 116.2581

## Evaluate I^2 to check if IVs appopirate for MR-Egger (Bowden et al.)

CRP_IVs_Isq <-Isq(R2_input$beta.exposure, 
                         R2_input$se.exposure)

## I^2 = 0.9942413

## Perform MR

Univariable_CRP <- mr(CRP_pneu_harm,
                      method_list = c("mr_ivw_mre",
                                      "mr_ivw_fe",
                                      "mr_weighted_median",
                                      "mr_weighted_mode",
                                      "mr_egger_regression"))

OR_CRP_to_pneumonia <- generate_odds_ratios(Univariable_CRP)

write.table(OR_CRP_to_pneumonia, file="Univariable_CRP_to_pneumonia.txt",
            sep = "\t", row.names = F, quote = F)

## Perform senstivity analyses

## Heterogeneity test via Cochran's Q

CRP_pneumonia_het <- mr_heterogeneity(CRP_pneu_harm)

## Test if MR-Egger intercept is significantly different from zero

CRP_pneumonia_Egger_intercept <- mr_pleiotropy_test(CRP_pneu_harm)

## Leave one out analysis using IVW-MRE and MR-Egger

CRP_pneumonia_LOO_IVW <- mr_leaveoneout(CRP_pneu_harm, method = mr_ivw_mre)

CRP_pneumonia_LOO_IVW %>% filter(p > 0.05)

CRP_pneumonia_LOO_Egger <- mr_leaveoneout(CRP_pneu_harm, method = mr_egger_regression)

CRP_pneumonia_LOO_Egger %>% filter(p > 0.05)

CRP_single_SNP <- mr_singlesnp(CRP_pneu_harm)

## Derive top three single SNPs

CRP_single_SNP <-  CRP_single_SNP[order(CRP_single_SNP$p),]

mr_forest_plot(CRP_single_SNP, exponentiate = F)


OR_CRP_to_pneumonia$METHOD <- as.factor(c("IVW-MRE", "IVW-FE", "Weighted Median", "Weighted Mode", "MR-Egger"))

## Forest plot of different MR methods effect size

FP_MR <- ggplot(data = OR_CRP_to_pneumonia, aes(x=METHOD, y=or, ymin=or_lci95, ymax=or_uci95)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  theme_bw() + 
  ylab("OR [95% CI] per SD increase in C-reactive protein") +
  xlab("MR model")


Filtered_CRP <- CRP_pneu_harm %>% filter(mr_keep == TRUE)

MR_PRESSO_input_CRP <- Filtered_CRP[,c("SNP", "beta.outcome", "beta.exposure", "se.outcome", "se.exposure")]

MR_PRESSO_CRP <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
                           DISTORTIONtest = TRUE, data = MR_PRESSO_input_CRP,
                           NbDistribution = 10000, SignifThreshold = 0.05)


