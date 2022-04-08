####################################

## Observational effect of GGT on pneumonia

## William Reay (2020)

####################################


library(dplyr)
library(data.table)
library(easyGgplot2)
library(car)
library(cowplot)
library(readr)
library(nnet)


## Load R data frame

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Pneumonia_phenotypes/")

source("Pneumonia_phenotypes.r")

## Import self-reported OR ICD10 pneumonia phenotype

Self_report_or_ICD10 <- fread("../Self_reported_OR_ICD10_pneumonia.tab.txt", header = T)

## Merged

UKBB_self_report <- merge(bd, Self_report_or_ICD10, by="f.eid")

## Import ICD-10 pneumonia phenotype

ICD_10 <- fread("../ICD10_pneumonia_self_reported_pneumonia_removed_as_controls.tab.txt", header = T)

UKBB_self_report_ICD10 <- merge(UKBB_self_report, ICD_10, by="f.eid")

UKBB_self_report_ICD10$SEX <- as.factor(ifelse(UKBB_self_report_ICD10$f.31.0.0 == "Male", 1, 0))

UKBB_self_report_ICD10$SMOKING <- as.factor(ifelse(UKBB_self_report_ICD10$f.20160.0.0=="Yes", 1, 0))

UKBB_self_report_ICD10$Self_reported_or_ICD <- as.factor(UKBB_self_report_ICD10$Self_reported_or_ICD)

UKBB_self_report_ICD10$ICD10_pneumonia <- as.factor(UKBB_self_report_ICD10$ICD10_pneumonia)


################ GGT estimate ###################

## Filter missing values for TG

GGT_UKBB_self_report <- UKBB_self_report_ICD10 %>% filter(!is.na(f.30730.0.0))

## Code individuals as in the top decile (GGT > 66.9 U/L)

GGT_UKBB_self_report$High_GGT <- ifelse(GGT_UKBB_self_report$f.30730.0.0 > 66.9, 1, 0)

## Scale TG to have zero mean and unit variance

GGT_UKBB_self_report$GGT_scaled <- as.numeric(scale(GGT_UKBB_self_report$f.30730.0.0))

GGT_UKBB_self_report$Pneumonia <- ifelse(GGT_UKBB_self_report$ICD10_pneumonia == "1", "Pneumonia", "Control")

## Test association - ICD10 pneumonia - TG as continous variable

ICD10_GGT_UKBB <- glm(ICD10_pneumonia ~ SEX*f.21003.0.0 + SEX*I(f.21003.0.0^2) +  SMOKING + f.189.0.0 + f.21003.0.0 + GGT_scaled, 
                     family="binomial", data = GGT_UKBB_self_report)


exp(cbind(coef(ICD10_GGT_UKBB), confint(ICD10_GGT_UKBB, level = 0.95)))

## Test association - ICD10 pneumonia - top decile GGT

Decile_ICD10_GGT_UKBB <- glm(ICD10_pneumonia ~ SEX*f.21003.0.0 + SEX*I(f.21003.0.0^2) +  SMOKING + f.189.0.0 + f.21003.0.0 + High_GGT, 
                            family="binomial", data = GGT_UKBB_self_report)

exp(cbind(coef(Decile_ICD10_GGT_UKBB), confint(Decile_ICD10_GGT_UKBB)))

ggplot2.boxplot(data=GGT_UKBB_self_report, xName='Pneumonia', yName='f.30730.0.0', groupName = 'Pneumonia',
                outlier.size=0.8, outlier.alpha = 0.2, ytitle="Gamma-glutamyltransferase (U/L)",
                xtitle="ICD-10 pneumonia diagnosis", groupColors=c('#999999','#E69F00'),
                showLegend = FALSE, backgroundColor="white", xtitleFont=c(12,"bold", "black"),
                ytitleFont=c(12,"bold", "black"), xTickLabelFont=c(10, "plain", "black"),
                yTickLabelFont=c(10, "plain", "black"), ylim=c(0, 250))

pvalue.extreme.z <- function(z) {
  log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
  mantissa <- 10^(log10.pvalue %% 1)
  exponent <- log10.pvalue %/% 1
  return(list(mantissa=mantissa,exponent=exponent))
}


## Non-smoking, women under 45 - triglycerides

TG_young_non_smoking_women <- TG_UKBB_self_report %>% filter(f.21003.0.0 <= 45 & SMOKING == "0" & SEX == "0")

Young_non_smokers_TG <- glm(ICD10_pneumonia ~ TG_scaled + f.189.0.0, family="binomial", data = TG_young_non_smoking_women)

exp(cbind(coef(Young_non_smokers_TG), confint(Young_non_smokers_TG)))

## Non-smoking, women under 45 - GGT

GGT_young_non_smoking_women <- GGT_UKBB_self_report %>% filter(f.21003.0.0 <= 45 & SMOKING == "0" & SEX == "0")

Young_non_smokers_GGT <- glm(ICD10_pneumonia ~ GGT_scaled + f.189.0.0, family="binomial", data = GGT_young_non_smoking_women)

exp(cbind(coef(Young_non_smokers_GGT), confint(Young_non_smokers_GGT)))



