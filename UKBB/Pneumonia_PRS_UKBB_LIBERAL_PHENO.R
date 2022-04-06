#############################################

## Pneumonia PRS - 80/20 split

## William Reay (2022)

## Clinically ascertained and self-reported

## MHC included and MHC excluded 

#############################################


## Define test and training set 

library(dplyr)
library(data.table)
library(rcompanion)
library(caret)
library(easyGgplot2)

set.seed(3636223)

Covariates <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Covariates_pneumonia_UKBB_white_British.tab.txt", header = T)

## Load broad pheno def

Broad_pneumonia <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Self_reported_or_ICD10_in_genetics.tab.txt", header = T)

## Merge

Merged_Broad <- merge(Covariates, Broad_pneumonia, by = "IID")

## Derive 80-20, training-testing split 

trainIndex <- createDataPartition(Merged_Broad$Self_reported_or_ICD, p = 0.8, list = FALSE)

Broad_train <- Merged_Broad[trainIndex, ]
Broad_test <- Merged_Broad[-trainIndex, ]

## Training 12024 cases and 256673 controls, Testing 3134 cases and 63940 controls

## Test MHC included first

Broad_train$Batch <- as.factor(Broad_train$Batch)
Broad_test$Batch <- as.factor(Broad_test$Batch)

MHC_incl <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/C_T_MHC_included.txt")

Colnames_MHC_incl <- make.unique(names(MHC_incl))

colnames(MHC_incl) <- Colnames_MHC_incl

MHC_incl <- MHC_incl %>% select(ends_with(c("IID", "AVG")))

MHC_incl_merged_TRAIN <- merge(MHC_incl, ICD10_train, by = "IID")

colnames(MHC_incl_merged_TRAIN)[2:9] <- paste("PRS", colnames(MHC_incl_merged_TRAIN[,c(2:9)]), sep = "_")

## Scale scores to have mean = zero and sd = 1

MHC_incl_merged_TRAIN[,c(2:9)] <- lapply(MHC_incl_merged_TRAIN[,c(2:9)], function(x) c(scale(x)))

## Define function to automate baseline testing

Scores_to_test <- as.list(colnames(MHC_incl_merged_TRAIN[,c(2:9)]))

PES_PRS_UKBB <- function(v, df){
  Score_model <- glm(glue::glue('Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                                PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + {v} + Batch'), family = "binomial", 
                     data = df)
  return(summary(Score_model)) 
}

Pneumonia_score <- sapply(Scores_to_test, PES_PRS_UKBB, df = MHC_incl_merged_TRAIN)

## Extract beta, se, z, and p value for each gene

Pneumonia_extract <- apply(Pneumonia_score, 2, function(x) return(as.data.frame(x$coefficients)[24,1:4]))

Pneumonia_results <- data.frame()

for (i in 1:length(Pneumonia_extract)) {
  Pneumonia_results <- rbind(Pneumonia_results, Pneumonia_extract[[i]])
}

Pneumonia_results$Score <- unlist(Scores_to_test)

write.csv(Pneumonia_results, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/PRS_results/BROAD_PHENO_MHC_included_PRS_TRAIN.csv",
          row.names = F, quote = F)

## Get variance on the liability scale (12%)

Pneumonia_PES_PRS_r2 <- function(v, k, p, df){
  Null_model <- glm(Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                      PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch, family = "binomial", data = df)
  Score_model <- glm(glue::glue('Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                                PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch + {v}'), family = "binomial", data = df)
  R2 <- nagelkerke(Score_model, null = Null_model)
  ## c.f. https://github.com/kn3in/ABC/blob/master/functions.R
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  R2_liab <- cc * e * R2$Pseudo.R.squared.for.model.vs.null[3, ] / (1 + cc * e* theta * R2$Pseudo.R.squared.for.model.vs.null[3, ])
  return(R2_liab)
}

Pneumonia_train_r2 <- sapply(Scores_to_test, Pneumonia_PES_PRS_r2, k=0.1296, p=0.19, df=MHC_incl_merged_TRAIN)

R2_liability_results <- data.frame()

for (i in 1:length(Pneumonia_train_r2)) {
  R2_liability_results <- rbind(R2_liability_results, Pneumonia_train_r2[[i]])
}

R2_liability_results$Score <- unlist(Scores_to_test)

write.csv(R2_liability_results, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/PRS_results/BROAD_PHENO_Pneumonia_MHC_incl_train_R2_conventional_weighting_liability.csv",
          row.names = F, quote = F)

## MHC excluded in training

MHC_excl <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/C_T_MHC_excluded.txt")

Colnames_MHC_excl <- make.unique(names(MHC_excl))

colnames(MHC_excl) <- Colnames_MHC_excl

MHC_excl <- MHC_excl %>% select(ends_with(c("IID", "AVG")))

MHC_excl_merged_TRAIN <- merge(MHC_excl, ICD10_train, by = "IID")

colnames(MHC_excl_merged_TRAIN)[2:9] <- paste("PRS", colnames(MHC_excl_merged_TRAIN[,c(2:9)]), sep = "_")

## Scale scores to have mean = zero and sd = 1

MHC_excl_merged_TRAIN[,c(2:9)] <- lapply(MHC_excl_merged_TRAIN[,c(2:9)], function(x) c(scale(x)))

Pneumonia_score_MHC_excl <- sapply(Scores_to_test, PES_PRS_UKBB, df = MHC_excl_merged_TRAIN)

## Extract beta, se, z, and p value for each gene

MHC_excl_Pneumonia_extract <- apply(Pneumonia_score_MHC_excl, 2, function(x) return(as.data.frame(x$coefficients)[24,1:4]))

MHC_excl_Pneumonia_results <- data.frame()

for (i in 1:length(MHC_excl_Pneumonia_extract)) {
  MHC_excl_Pneumonia_results <- rbind(MHC_excl_Pneumonia_results, MHC_excl_Pneumonia_extract[[i]])
}

MHC_excl_Pneumonia_results$Score <- unlist(Scores_to_test)

write.csv(MHC_excl_Pneumonia_results, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/PRS_results/BROAD_PHENO_MHC_excluded_PRS_TRAIN.csv",
          row.names = F, quote = F)

MHC_excl_Pneumonia_train_r2 <- sapply(Scores_to_test, Pneumonia_PES_PRS_r2, k=0.1296, p=0.19, df=MHC_excl_merged_TRAIN)

MHC_excl_R2_liability_results <- data.frame()

for (i in 1:length(MHC_excl_Pneumonia_train_r2)) {
  MHC_excl_R2_liability_results <- rbind(MHC_excl_R2_liability_results, MHC_excl_Pneumonia_train_r2[[i]])
}

MHC_excl_R2_liability_results$Score <- unlist(Scores_to_test)

write.csv(MHC_excl_R2_liability_results, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/UKBB/PRS_results/BROAD_Pheno_Pneumonia_MHC_excl_train_R2_conventional_weighting_liability.csv",
          row.names = F, quote = F)

## Evaluate best performing score (Pt < 0.1) in test set and get OR for each decile relative to 5-6

Merged_TEST <- merge(MHC_incl, Broad_test, by = "IID")

Merged_TEST$Test_PRS <- Merged_TEST$`0_1_AVG`

Merged_TEST$Test_PRS <- as.numeric(scale(Merged_TEST$Test_PRS))

Test_effect <- glm(Self_reported_or_ICD ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                     PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch + Test_PRS,
                   family = "binomial", data = Merged_TEST)


Merged_TEST$Decile <- as.factor(ntile(Merged_TEST$Test_PRS, 10))

Decile_test <- Test_effect <- glm(Self_reported_or_ICD~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 +
                                    PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + Batch + Decile,
                                  family = "binomial", data = Merged_TEST)


Broad_dens <- ggplot2.density(data=Merged_TEST, 
                              xName='Test_PRS', groupName='Self_reported_or_ICD', 
                              alpha=0.5, fillGroupDensity=TRUE, addMeanLine=FALSE, meanLineColor="black", xtitle="Pneumonia PRS (SD)", 
                              ytitle="Density", xtitleFont=c(12, "plain", "black"), ytitleFont=c(12, "plain", "black"), 
                              axisLine=c(0.5, "solid", "black"), xTickLabelFont=c(12, "plain", "black"), 
                              yTickLabelFont=c(12, "plain", "black"), backgroundColor="white", removePanelGrid=TRUE, faceting=TRUE, bins=30, 
                              densityAlpha=0.3, densityLineSize=4, groupColors=c('#999999','#0099FF'), mainTitle="Clinically ascertained or self-reported pneumonia",
                              legendPosition="none")
