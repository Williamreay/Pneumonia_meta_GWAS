## MR-pheWAS script - IEUGWAS

## William Reay (2021)


library(data.table)
library(ieugwasr)
library(dplyr)

## Picked the gene from each region with an eQTL with the highest PIP for the credible set from the eQTL catalog


IEU_mr_pheWAS <- function(SNP, EA, NEA, IV_beta) {
  
  ## Extract IV-outcome effect sizes from the IEUGWAS db
  
  SNP_pheWAS <- phewas(SNP, pval = 1)
  
  ## Check for mismatching effect and non-effect alleles
  
  mismatch_SNP <- which(SNP_pheWAS$ea != EA & SNP_pheWAS$nea != NEA,arr.ind=TRUE)
  SNP_pheWAS[mismatch_SNP,]$beta <- SNP_pheWAS[mismatch_SNP,]$beta*-1
  SNP_pheWAS[mismatch_SNP,]$ea <- EA
  SNP_pheWAS[mismatch_SNP,]$nea <- NEA
  
  ## Perform MR
  
  SNP_pheWAS$MR_beta <- SNP_pheWAS$beta/IV_beta
  SNP_pheWAS$MR_SE <- abs(SNP_pheWAS$se/IV_beta)
  SNP_pheWAS$MR_pval <- pnorm(abs(SNP_pheWAS$MR_beta)/SNP_pheWAS$MR_SE, lower.tail=FALSE) * 2
  
  ## Retain the following - GWAS catalog imports, Brain regions, IEUGWAS UKBB analyses, Shin et al metabolites, Kettunen metabolites, UKBB metabolites, immune metabolites, and curated consortia
  
  SNP_pheWAS <- dplyr::filter(SNP_pheWAS, grepl('ukb-b|met-a|met-b|met-c|met-d|ubm-a|ebi-a|ieu-a|ieu-b', id))

  return(SNP_pheWAS)
}

## TNFRSF1A IV - rs1800692, G allele, beta = -0.071, SE = 0.0088, PIP = 0.98 (Lepik 2017 blood)

TNFRSF1A_MR <- IEU_mr_pheWAS("rs1800692", "G", "A", -0.071)

## Remove incorrectly formatted sumstats from ieugwas

TNFRSF1A_MR <- TNFRSF1A_MR %>% filter(id != "ebi-a-GCST007799" & id != "ebi-a-GCST007800")

TNFRSF1A_MR$MR_FWER <- p.adjust(TNFRSF1A_MR$MR_pval, method="bonferroni")
TNFRSF1A_MR$MR_FDR <- p.adjust(TNFRSF1A_MR$MR_pval, method="fdr")

write.table(TNFRSF1A_MR, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MR_pheWAS_drug_targets/TNFRSF1A_ieugwas_MR_PheWAS.txt",
            sep = "\t", row.names = F, quote = F)

## MUC5AC IV - C allele, beta = -0.218184, rsID = 	rs28515631

MUC5AC_MR <- IEU_mr_pheWAS("rs28515631", "C", "A", -0.218184)

MUC5AC_MR <- MUC5AC_MR %>% filter(id != "ebi-a-GCST007799" & id != "ebi-a-GCST007800")

MUC5AC_MR$MR_FWER <- p.adjust(MUC5AC_MR$MR_pval, method="bonferroni")
MUC5AC_MR$MR_FDR <- p.adjust(MUC5AC_MR$MR_pval, method="fdr")

write.table(MUC5AC_MR, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MR_pheWAS_drug_targets/MUC5AC_ieugwas_MR_PheWAS.txt",
            sep = "\t", row.names = F, quote = F)
