#####################################

## Pneumonia heterogeneity vs significance

## William Reay (2022)

#####################################

library(data.table)
library(ggplot2)
library(dplyr)

Pneu <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/CHR_BP_annotated_FINAL_IVW_meta.txt.gz", header = T)

Pneu <- rename(Pneu,"SNP"="rsID")


Pneu$logGWASP <- -log10(Pneu$P)
Pneu$logHetP <- -log10(Pneu$HetPVal)


HapMap <- fread("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/Munged/IVW_META_munged_HapMap3.sumstats.gz", header = T)

Pneu <- merge(Pneu, HapMap, by ="SNP")

Pneu <- Pneu %>% filter(P < 0.05)

ggplot(data = Pneu,
       aes(x = HetISq, y=-log10(P))) +
  geom_point(alpha = 0.5, size = 1.1) +
  geom_vline(xintercept = 75, linetype="longdash", colour="blue") +
  geom_hline(yintercept = 5, linetype="longdash", colour="blue") +
  xlab("I-squared statistic") +
  ylab("-log10 Pneumonia P-value") +
  ggtitle("I-squared statistic vs GWAS P-value")

ggplot(data = Pneu, aes(x = -log10(HetPVal), y=-log10(P))) +
  geom_point(alpha = 0.5, size = 1.1) +
  geom_vline(xintercept = 1.3, linetype="longdash", colour="red") +
  geom_hline(yintercept = 5, linetype="longdash", colour="red") +
  theme_bw() +
  xlab("-log10 Cochran's Q P-value") +
  ylab("-log10 Pneumonia P-value") +
  ggtitle("Heterogeneity P-value vs GWAS P-value")
