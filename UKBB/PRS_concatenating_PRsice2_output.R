#################################

## Concatenating chromosome-wide pneumonia PRS at each P-value threshold

## William Reay (2020)

#################################

library(data.table)
library(dplyr)
library(plyr)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/")

## Import PRS at each P threshold on all 22 autosomes, concatenated into the sample file

All_chr_PRS <- fread("PRS_RAW/PRS_all_chr_all_score.txt", sep = " ", header = T)

Colnames_all_chr_PRS <- make.unique(names(All_chr_PRS))

colnames(All_chr_PRS) <- Colnames_all_chr_PRS


## Convert double IID to single IID type
  
All_chr_PRS$IID <- gsub(".*_","",All_chr_PRS$IID)  

## Select each P value threshold
  
## P < 1e-05
  
P_0.00001_threshold <- All_chr_PRS %>% select(IID, starts_with("Pt_1e-05"))

## P < 0.005

P_0.005_threshold <- All_chr_PRS %>% select(IID, starts_with("Pt_0.005"))

## P < 0.05

P_0.05_threshold <- All_chr_PRS %>% select(IID, starts_with("Pt_0.05"))

## P < 0.5

P_0.5_threshold <- All_chr_PRS %>% select(IID, starts_with("Pt_0.5"))

## P < 1

P_1_threshold <- All_chr_PRS %>% select(IID, starts_with("Pt_1."), ends_with("Pt_1"))


## Sum each chrosome to obtain a genome-wide PRS at the specified threshold

## P < 1e-05

P_0.00001_threshold$PRS_0.00001 <- rowSums(P_0.00001_threshold[,-1])

P_0.00001_threshold <- P_0.00001_threshold %>% select(IID, PRS_0.00001)

## P < 0.005

P_0.005_threshold$PRS_0.005 <- rowSums(P_0.005_threshold[,-1])

P_0.005_threshold <- P_0.005_threshold %>% select(IID, PRS_0.005)

## P < 0.05

P_0.05_threshold$PRS_0.05 <- rowSums(P_0.05_threshold[,-1])

P_0.05_threshold <- P_0.05_threshold %>% select(IID, PRS_0.05)

## P < 0.5

P_0.5_threshold$PRS_0.5 <- rowSums(P_0.5_threshold[,-1])

P_0.5_threshold <- P_0.5_threshold %>% select(IID, PRS_0.5)

## P < 1

P_1_threshold$PRS_1 <- rowSums(P_1_threshold[,-1])

P_1_threshold <- P_1_threshold %>% select(IID, PRS_1)

## Merge all PRS dfs together

All_PRS <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"),
                           list(P_0.00001_threshold, P_0.005_threshold, P_0.05_threshold, P_0.5_threshold, P_1_threshold))

## Output

write.table(All_PRS, file="PRS_scores/All_RAW_PRS_scores_all_chromosome.tab.txt",
            sep = "\t", row.names = F, quote = F)

