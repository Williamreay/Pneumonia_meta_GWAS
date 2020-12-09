################################

## Preparing antigen data

## William Reay (2020)

#################################

library(dplyr)
library(DescTools)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/")


## Read in seropositivity data


RAW_antigen_dat <- fread("Antigen_data/UKBB_seropositivity.txt", header = T)

NOMISS_antigen_dat <- RAW_antigen_dat %>% filter(!is.na(f.23050.0.0) | !is.na(f.23051.0.0) | 
                                                   !is.na(f.23052.0.0) | !is.na(f.23053.0.0) | 
                                                   !is.na(f.23054.0.0) | !is.na(f.23055.0.0) | 
                                                   !is.na(f.23056.0.0) | !is.na(f.23057.0.0) | 
                                                   !is.na(f.23058.0.0) | !is.na(f.23059.0.0) | 
                                                   !is.na(f.23060.0.0) | !is.na(f.23061.0.0) | 
                                                   !is.na(f.23062.0.0) | !is.na(f.23063.0.0) | 
                                                   !is.na(f.23064.0.0) | !is.na(f.23065.0.0) | 
                                                   !is.na(f.23066.0.0) | !is.na(f.23067.0.0) | 
                                                   !is.na(f.23068.0.0) | !is.na(f.23069.0.0) | 
                                                   !is.na(f.23070.0.0) | !is.na(f.23071.0.0) | 
                                                   !is.na(f.23073.0.0) | !is.na(f.23074.0.0))



## Load covariate data to merge with IDs with genetics available

NOMISS_antigen_dat <- rename(NOMISS_antigen_dat, "IID"="f.eid")

Covariates <- fread("Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)                                                            

Gen_IDs_merged <- merge(NOMISS_antigen_dat, Covariates, by="IID")

## Load PRS and PES

Complement_PES <- fread("PES_scores/Complement_pathway/Complement_pathway_pneumonia_PES.tab.txt", header = T)
Bile_PES <- fread("PES_scores/Bile_acid_metabolism/Bile_acid_metabolism_pneumonia_PES.tab.txt", header = T)
Lectin_PES <- fread("PES_scores/Lectin_pathway/Lectin_pathway_pneumonia_PES.tab.txt", header = T)
p53_PES <- fread("PES_scores/p53_pathway/p53_pathway_pneumonia_PES.tab.txt", header = T)
RIG1_PES <- fread("PES_scores/RIG_1_pathway/RIG1_pathway_pneumonia_PES.tab.txt", header = T)

PRS_all <- fread("PRS_scores/All_RAW_PRS_scores_all_chromosome.tab.txt")

PRS_all <- PRS_all %>% select(IID, starts_with("PRS"))

## Merge the PRS and PES into the same df

Merged_PES_PRS_antigen <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                                 list(PRS_all, Complement_PES, Bile_PES, Lectin_PES, p53_PES, RIG1_PES, Gen_IDs_merged))

## Make list of columes to apply function over

MFI_value_col <- as.list(colnames(Merged_PES_PRS_antigen %>% select(starts_with("f.2305") | starts_with("f.2306") |
                                                                      starts_with("f.2307"))))

PES_PRS_colnames <- as.list(colnames(Merged_PES_PRS_antigen %>% select(starts_with("PRS") | ends_with("PES"))))


## Define QC covariates - potential spillover and one additional freeze and thaw cycle

Merged_PES_PRS_antigen$Spill_over <- ifelse(is.na(Merged_PES_PRS_antigen$f.23049.0.0) | Merged_PES_PRS_antigen$f.23049.0.0  == 2, 0, 1)

Merged_PES_PRS_antigen$Freeze_thaw <- ifelse(is.na(Merged_PES_PRS_antigen$f.23049.0.0) | Merged_PES_PRS_antigen$f.23049.0.0  == 1, 0, 1)

## Save df as RDS

saveRDS(Merged_PES_PRS_antigen, file="Antigen_data/Cleaned_antigen_data.rds")

## Make df with triglycerides and GGT

source("../../Pneumonia_phenotypes/Pneumonia_phenotypes.r")


## Filter missing values for TG or GGT

TG_UKBB <- bd %>% filter(!is.na(f.30870.0.0))
TG_UKBB <- TG_UKBB %>% filter(!is.na(f.30730.0.0))

## Scale TG and GGT

TG_UKBB$Scaled_TG <- as.numeric(scale(TG_UKBB$f.30870.0.0))
TG_UKBB$Scaled_GGT <- as.numeric(scale(TG_UKBB$f.30730.0.0))

## Winsorize GGT and TG

TG_UKBB$Unscaled_winsorized_TG <- Winsorize(TG_UKBB$f.30870.0.0, min = 0, max = 4.829708)
TG_UKBB$Unscaled_winsorized_GGT <- Winsorize(TG_UKBB$f.30870.0.0, min = 0, max = 163.3996)

TG_UKBB$Scaled_winsorized_TG <- as.numeric(scale(TG_UKBB$Unscaled_winsorized_TG))
TG_UKBB$Scaled_winsorized_GGT <- as.numeric(scale(TG_UKBB$Unscaled_winsorized_GGT))

GGT_TG_output <- TG_UKBB %>% select(f.eid, Scaled_GGT, Scaled_TG, Scaled_winsorized_GGT, Scaled_winsorized_TG)

GGT_TG_output <- rename(GGT_TG_output, "IID"="f.eid")

Antigen_GGT_TG <- merge(Merged_PES_PRS_antigen, GGT_TG_output, by= "IID")

saveRDS(Antigen_GGT_TG, file="Antigen_data/Antigen_GGT_TG.rds")
