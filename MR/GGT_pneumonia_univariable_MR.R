########################################

## GGT ==> pneumonia susceptibility (2022)

## Univariable MR - total effect

## William Reay (2022)

########################################

library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(MRPRESSO)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/LCV/MR_GGT_follow_up/")

## Read in GGT exposure data in SD units

GGT_raw <- fread("CAUSE_formatted_GGT_30730_irnt.gwas.imputed_v3.both_sexes.tsv.gz")

## Retain GW sig var

GGT_gw_sig <- GGT_raw %>% filter(pval < 5e-08)

GGT_gw_sig$Phenotype <- "GGT"

## Format as IEUGWAS IV

GGT_exp <- format_data(dat = GGT_gw_sig,
                      type = "exposure",
                      beta_col = "beta",
                      se_col = "se",
                      snp_col = "rsid",
                      effect_allele_col = "alt",
                      other_allele_col = "ref",
                      phenotype_col = "Phenotype")

## LD clump per usual r2

GGT_clumped <- clump_data(GGT_exp, pop = "EUR")

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

GGT_pneu_harm <- harmonise_data(GGT_clumped, Pneumonia_outcome, action = 3)


## Estimate variance explained and calculate i) F-stat, and ii) I^2 stat

R2_input <- GGT_pneu_harm %>% filter(mr_keep == TRUE)

Approximated_r2_gamma_glutamyltransferase <- sum(get_r_from_pn(R2_input$pval.exposure,344104)^2)

## Calculate F statistic to evaluate whether it's > 10

F_statistic <- function(r,n,k){r*(n-1-k)/(1-r)/k}

F_statistic(Approximated_r2_gamma_glutamyltransferase, 344104, 206)

## F stat = 114.7967

## Evaluate I^2 to check if IVs appopirate for MR-Egger (Bowden et al.)

Isq = function(y,s){
  k = length(y)
  w = 1/s^2; sum.w = sum(w)
  mu.hat = sum(y*w)/sum.w
  Q = sum(w*(y-mu.hat)^2)
  Isq = (Q - (k-1))/Q
  Isq = max(0,Isq)
  return(Isq)
}

Gamma_glut_IVs_Isq <-Isq(R2_input$beta.exposure, 
                         R2_input$se.exposure)

## I^2 = 0.991847

## Perform MR

Univariable_GGT <- mr(GGT_pneu_harm,
                      method_list = c("mr_ivw_mre",
                                      "mr_ivw_fe",
                                      "mr_weighted_median",
                                      "mr_weighted_mode",
                                      "mr_egger_regression"))

OR_GGT_to_pneumonia <- generate_odds_ratios(Univariable_GGT)

write.table(OR_GGT_to_pneumonia, file="Univariable_GGT_to_pneumonia.txt",
            sep = "\t", row.names = F, quote = F)

## Perform senstivity analyses

## Heterogeneity test via Cochran's Q

GGT_pneumonia_het <- mr_heterogeneity(GGT_pneu_harm)

## Test if MR-Egger intercept is significantly different from zero

GGT_pneumonia_Egger_intercept <- mr_pleiotropy_test(GGT_pneu_harm)

## Leave one out analysis using IVW-MRE and MR-Egger

GGT_pneumonia_LOO_IVW <- mr_leaveoneout(GGT_pneu_harm, method = mr_ivw_mre)

GGT_pneumonia_LOO_IVW %>% filter(p < 0.1)

GGT_pneumonia_LOO_Egger <- mr_leaveoneout(GGT_pneu_harm, method = mr_egger_regression)

GGT_pneumonia_LOO_Egger %>% filter(p < 0.1)

OR_GGT_to_pneumonia$METHOD <- as.factor(c("IVW-MRE", "IVW-FE", "Weighted Median", "Weighted Mode", "MR-Egger"))

## Forest plot of different MR methods effect size

FP_MR <- ggplot(data = OR_GGT_to_pneumonia, aes(x=METHOD, y=or, ymin=or_lci95, ymax=or_uci95)) +
  geom_pointrange() +
  geom_hline(yintercept=1, lty=2) +
  coord_flip() +
  theme_bw() + 
  ylab("OR [95% CI] per SD increase in gamma glutamyltransferase") +
  xlab("MR model")


Filtered_GGT <- GGT_pneu_harm %>% filter(mr_keep == TRUE)

MR_PRESSO_input_GGT <- Filtered_GGT[,c("SNP", "beta.outcome", "beta.exposure", "se.outcome", "se.exposure")]

MR_PRESSO_GGT <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE,
                           DISTORTIONtest = TRUE, data = MR_PRESSO_input_GGT,
                           NbDistribution = 10000, SignifThreshold = 0.05)


