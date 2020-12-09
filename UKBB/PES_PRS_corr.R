###############################################

## PES ~ PRS correlation in UKBB

## William Reay (2020)

###############################################


library(corrplot)
library(dplyr)
library(data.table)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/")

## Load raw, unscaled, PRS

PRS_all <- fread("PRS_scores/All_RAW_PRS_scores_all_chromosome.tab.txt")

PRS_all <- PRS_all %>% select(IID, starts_with("PRS"))

## Load in PES

Complement_PES <- fread("PES_scores/Complement_pathway/Complement_pathway_pneumonia_PES.tab.txt", header = T)
Bile_PES <- fread("PES_scores/Bile_acid_metabolism/Bile_acid_metabolism_pneumonia_PES.tab.txt", header = T)
Lectin_PES <- fread("PES_scores/Lectin_pathway/Lectin_pathway_pneumonia_PES.tab.txt", header = T)
p53_PES <- fread("PES_scores/p53_pathway/p53_pathway_pneumonia_PES.tab.txt", header = T)
RIG1_PES <- fread("PES_scores/RIG_1_pathway/RIG1_pathway_pneumonia_PES.tab.txt", header = T)

IID <- Bile_PES[, 1]
Complement_PES <- Complement_PES[, 8]
Bile_PES <- Bile_PES[, 2]
Lectin_PES <- Lectin_PES[, 7]
p53_PES <- p53_PES[, 20]
RIG1_PES <- RIG1_PES[, 9]


PES_combined <- cbind(IID, Complement_PES, Bile_PES, Lectin_PES, p53_PES, RIG1_PES)

## Merge PES and PRS

PES_PRS <- merge(PRS_all, PES_combined, by="IID")

PES_PRS$PRS_0.00001 <- as.numeric(scale(PES_PRS$PRS_0.00001))
PES_PRS$PRS_0.005 <- as.numeric(scale(PES_PRS$PRS_0.005))
PES_PRS$PRS_0.05 <- as.numeric(scale(PES_PRS$PRS_0.05))
PES_PRS$PRS_0.5 <- as.numeric(scale(PES_PRS$PRS_0.5))
PES_PRS$PRS_1 <- as.numeric(scale(PES_PRS$PRS_1))
PES_PRS$Bile_acid_PES <- as.numeric(scale(PES_PRS$Bile_acid_PES))
PES_PRS$Complement_pathway_PES <- as.numeric(scale(PES_PRS$Complement_pathway_PES))
PES_PRS$Lectin_pathway_PES <- as.numeric(scale(PES_PRS$Lectin_pathway_PES))
PES_PRS$p53_pathway_PES <- as.numeric(scale(PES_PRS$p53_pathway_PES))
PES_PRS$RIG1_pathway_PES <- as.numeric(scale(PES_PRS$RIG1_pathway_PES))

write.table(PES_PRS, file="Scaled_PES_PRS_UKBB_pneumonia.tab.txt",
            sep="\t", row.names = F, quote = F)

Corr_PES_PRS <- cor(PES_PRS[, -1])


col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(Corr_PES_PRS, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.45, number.cex= 0.5)

## Generate correlation plot which only shows significant correlations

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(Corr_PES_PRS)

corrplot(Corr_PES_PRS, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.0005, insig = "blank")

corrplot(Corr_PES_PRS, method = "color", addCoef.col = "black", col = col(200),
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.6, number.cex= 0.5, p.mat = p.mat, sig.level = 0.05, insig = "blank")
