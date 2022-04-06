########################################

## MR-Clust - CRP ==> Pneumonia

## William Reay (2022)


########################################


library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)
library(mrclust)

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

## Delineate bx and by effect size/SE and ratio estimates for MRClust input

bx = CRP_pneu_harm$beta.exposure
bxse =  CRP_pneu_harm$se.exposure

by = CRP_pneu_harm$beta.outcome
byse = CRP_pneu_harm$se.outcome

ratio_est = by/bx
ratio_est_se = byse/abs(bx)

snp_names = CRP_pneu_harm$SNP

## Default MR Clust

set.seed(2000030885)
res_em = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = snp_names)


names(res_em$results)

head(res_em$results$all)

head(res_em$results$best)

plot_sbp_best = res_em$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with CRP") +
  ggplot2::ylab("Genetic association with Pneumonia susceptibility") +
  ggplot2::ggtitle("")

clusters = unique(res_em$results$best$cluster_class)

## Variants located to clusters with allocation prob >0.8

res80 = mrclust::pr_clust(dta = res_em$results$best, prob = 0.8)


keep80 = which(snp_names %in% res80$observation)
bx80   = bx[keep80]
bxse80 = bxse[keep80]
by80   = by[keep80]
byse80 = byse[keep80]
snp_names80 = snp_names[keep80]
plot.sbp.pr80 = two_stage_plot(res = res80, bx = bx80, by = by80, bxse = bxse80,
                               byse = byse80, obs_names = snp_names80) + 
  ggplot2::xlim(0, max(abs(bx80) + 2*bxse80)) + 
  ggplot2::xlab("Genetic association with CRP") + 
  ggplot2::ylab("Genetic association with Pneumonia susceptibility") + 
  ggplot2::ggtitle("");
# plot result
plot.sbp.pr80



