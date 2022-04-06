###########################

## UKBB pneumonia covariate analyses

###########################

library(dplyr)
library(data.table)
library(ggplot2)
library(car)


## Load R data frame

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/Pneumonia_phenotypes/")

source("Pneumonia_phenotypes.r")

## Remove individuals who are missing the following:
#' Assessment age
#' Sex
#' Townsend Deprivation index
#' Smoking status

Missing_cov_removed <- bd %>% filter(!is.na(f.31.0.0)) %>% filter(!is.na(f.21003.0.0)) %>% 
  filter(!is.na(f.189.0.0)) %>% filter(!is.na(f.20160.0.0))


### SELF-REPORTED PNEUMONIA ####

## Test association of base pneumonia covariates

## Define list of data columns

Self_report_cols <- as.list(paste("f.20002.",0:3, ".", 1:33,sep=""))

## Dummy code individuals with self-reported pneumonia (1398 coding) in any of the f.20002 cols == 1398

select_eid <- Missing_cov_removed %>%
  select(starts_with("f.eid"))

select_20002 <-  Missing_cov_removed %>%
  select(starts_with("f.20002."))

Self_reported_df <- bind_cols(select_eid, select_20002)

Self_reported_pneumonia <- which(Self_reported_df=="1398", arr.ind=T)

Self_reported_pneumonia_row <- Self_reported_pneumonia[,1]

Self_reported_df$Self_reported_pneumonia_col <- 0

Self_reported_df[Self_reported_pneumonia_row, ]$Self_reported_pneumonia_col <- 1

Self_reported_df$Self_reported_pneumonia_col <- as.factor(Self_reported_df$Self_reported_pneumonia_col)

Self_reported_df <- Self_reported_df %>% select(f.eid, Self_reported_pneumonia_col)

## Merge with original df

Pneumonia_self_report <- merge(Missing_cov_removed, Self_reported_df, by="f.eid")

## Test association with the following variables
#' Age
#' Sex
#' Townsend Deprivation Index
#' Smoking
#' 

Pneumonia_cov_test$SEX <- as.factor(ifelse(Pneumonia_cov_test$f.31.0.0 == "Male", 1, 0))

Pneumonia_cov_test$SMOKING <- as.factor(ifelse(Pneumonia_cov_test$f.20160.0.0=="Yes", 1, 0))

Pneumonia_cov <- glm(Self_reported_pneumonia_col ~ SMOKING + f.189.0.0 + SEX + f.21003.0.0, family = "binomial", data=Pneumonia_cov_test)


### ICD-10 codes ###

select_ICD10_primary <- Missing_cov_removed %>% select(starts_with("f.41202"))

select_ICD10_secondary <- Missing_cov_removed %>% select(starts_with("f.41204"))

ICD_df <- cbind(select_eid, select_ICD10_primary, select_ICD10_secondary)


## Derive list of pneumonia ICD codes for primary or secondary diagnosis 

## J100

J100_ICD <- which(ICD_df == "J100", arr.ind = T)

J100_ICD_row <- J100_ICD[, 1]

## J11

J11_ICD <- which(ICD_df == "J110", arr.ind = T)

J11_ICD_row <- J11_ICD[, 1]


## J12 

J12_ICD <- which(ICD_df == "J121" | ICD_df == "J122" | ICD_df == "J123" | ICD_df == "J128" | ICD_df == "J129", arr.ind = T)

J12_ICD_row <- J12_ICD[, 1]

## J13

J13_ICD <- which(ICD_df == "J13", arr.ind = T)

J13_ICD_row <- J13_ICD[, 1]

## J14

J14_ICD <- which(ICD_df == "J14", arr.ind = T)

J14_ICD_row <- J14_ICD[, 1]

## J15

J15_ICD <- which(ICD_df == "J150" | ICD_df == "J151" | ICD_df == "J152" | ICD_df == "J153" | ICD_df == "J154" | ICD_df == "J155" | ICD_df == "J156" | ICD_df == "J157" | ICD_df == "J158" | ICD_df == "J159", arr.ind = T)        

J15_ICD_row <- J15_ICD[, 1]

## J16

J16_ICD <- which(ICD_df == "J160" | ICD_df =="J168", arr.ind = T)

J16_ICD_row <- J16_ICD[, 1]

## J17

J17_ICD <- which(ICD_df == "J172" | ICD_df =="J173", arr.ind = T)

J17_ICD_row <- J17_ICD[, 1]

## J18

J18_ICD <- which(ICD_df == "J180" | ICD_df =="J181" | ICD_df =="J182" | ICD_df =="J188" | ICD_df =="J189", arr.ind = T)

J18_ICD_row <- J18_ICD[, 1]

## Pneumonia_ICD_total_rows

Pneumonia_total <- c(J100_ICD_row, J11_ICD_row, J12_ICD_row, J13_ICD_row, J14_ICD_row, J15_ICD_row, J16_ICD_row, J17_ICD_row, J18_ICD_row)

Pneumonia_total_uniq <- unique(Pneumonia_total)

## Code ICD-10 pneumonia column

ICD_df$ICD10_pneumonia <- 0

ICD_df[Pneumonia_total_uniq, ]$ICD10_pneumonia <- 1

ICD_df$ICD10_pneumonia <- as.factor(ICD_df$ICD10_pneumonia)

ICD10_df <- ICD_df %>% select(f.eid, ICD10_pneumonia)


## Merge with the non-missing covariate individuals

ICD_self_report_merged <- merge(ICD10_df, Pneumonia_self_report, by="f.eid")

ICD_merged$SEX <- as.factor(ifelse(ICD_merged$f.31.0.0 == "Male", 1, 0))

ICD_merged$SMOKING <- as.factor(ifelse(ICD_merged$f.20160.0.0=="Yes", 1, 0))

Pneumonia_cov_ICD10 <- glm(ICD10_pneumonia ~ SMOKING + f.189.0.0 + SEX + f.21003.0.0, family = "binomial", data=ICD_merged)

exp(cbind(coef(Pneumonia_cov_ICD10), confint(Pneumonia_cov_ICD10)))


## Define phenotypes

## ICD-10, all non-overlapping self-reported removed, i.e. individuals with self-reported pneumonia but without ICD10 codes

ICD_self_report_merged$Remove_controls <- ifelse(ICD_self_report_merged$ICD10_pneumonia == "0" & ICD_self_report_merged$Self_reported_pneumonia_col == "1", 1, 0)

ICD_self_report_merged_controls_removed <- ICD_self_report_merged %>% filter(Remove_controls == 0)

ICD_self_report_merged_controls_removed <- ICD_self_report_merged_controls_removed %>% select(f.eid, ICD10_pneumonia)

write.table(ICD_self_report_merged_controls_removed, file="../ICD10_pneumonia_self_reported_pneumonia_removed_as_controls.tab.txt",
            sep="\t", row.names = F, quote = F)

## ICD-10 OR Self-reported 

ICD_self_report_merged$Self_reported_or_ICD <- ifelse(ICD_self_report_merged$ICD10_pneumonia == "1" | ICD_self_report_merged$Self_reported_pneumonia_col == "1", 1, 0)

ICD_self_report_merged$Self_reported_or_ICD <- as.factor(ICD_self_report_merged$Self_reported_or_ICD)

ICD_self_report_merged <- ICD_self_report_merged %>% select(f.eid, Self_reported_or_ICD)

write.table(ICD_self_report_merged, file="../Self_reported_OR_ICD10_pneumonia.tab.txt",
            sep="\t", row.names = F, quote = F)
