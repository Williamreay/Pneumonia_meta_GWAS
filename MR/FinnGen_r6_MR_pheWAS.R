#################################

## FinnGen r6 - MR PheWAS

## William Reay (2021)

#################################

library(dplyr)
library(data.table)
library(ggplot2)

## Define function that includes derivation of SE from P val as not provided by FinnGen output
## Note that FinnGen already aligned to the ALT allele like the eQTL beta

FinnGen_MR_PheWAS <- function(SNP, EA, NEA, IV_beta, df) {
  
  df <- df %>% filter(beta != "NA")
  
  df$Z <- sign(df$beta) * sqrt(qchisq(df$pval, 1, lower =F))
  df$SE <- abs(df$beta/df$Z)
  
  ## MR - Wald ratio
  
  df$MR_beta <- df$beta/IV_beta
  df$MR_SE <- abs(df$SE/IV_beta)
  df$MR_pval <- pnorm(abs(df$MR_beta)/df$MR_SE, lower.tail=FALSE) * 2
  
  ## FWER and FDR correction
  
  df$FWER <- p.adjust(df$MR_pval, method="bonferroni")
  df$FDR <- p.adjust(df$MR_pval, method="fdr")
  
  return(df)
}

TNFRSF1A_df <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MR_pheWAS_drug_targets/12_6333180_A_G_phenotype_associations.tsv")

TNFRSF1A_MR_Finngen <- FinnGen_MR_PheWAS("rs1800692", "G", "A", -0.071, df = TNFRSF1A_df)

write.table(TNFRSF1A_MR_Finngen, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MR_pheWAS_drug_targets/TNFRSF1A_finngen_MR_PheWAS.txt",
            sep = "\t", row.names = F, quote = F)

