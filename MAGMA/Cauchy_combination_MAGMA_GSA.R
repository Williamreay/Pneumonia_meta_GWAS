###########################

## MAGMA GSA - Cauchy Combination

## William Reay (2022)

###########################

## Cauchy combination across all three

source("../Scripts/ACAT_function.R")

GSA_no <- fread("../Common_var_results/MAGMA/Multiple_testing_corr/gProfiler_sets_GSA/NO_bound_pneumonia_gprofiler.gsa.out", skip=3)
GSA_no <- GSA_no %>% filter(NGENES > 5)


GSA_cons <- fread("../Common_var_results/MAGMA/Multiple_testing_corr/gProfiler_sets_GSA/CONS_pneumonia_gprofiler.gsa.out", skip=3)
GSA_cons <- GSA_cons %>% filter(NGENES > 5)

GSA_lib <- fread("../Common_var_results/MAGMA/Multiple_testing_corr/gProfiler_sets_GSA/LIB_pneumonia_gprofiler.gsa.out", skip=3)
GSA_lib <- GSA_lib %>% filter(NGENES > 5)

Merged_GSA <- rbind(GSA_cons, GSA_lib, GSA_no)

genes <- sort(unique(Merged_GSA$VARIABLE))

gene_ps <- list()

for (i in genes) { gene_ps[[i]] <- c(Merged_GSA[Merged_GSA$VARIABLE==i,][,7])}

genes_acat <- list()
for (i in genes) { genes_acat[[i]] <- ACAT(unlist(gene_ps[[i]], use.names=FALSE))}

af_acat <- data.frame(cbind(as.vector(genes), as.vector(unlist(genes_acat))))

af_acat <- rename(af_acat, c("Gene_set_ID"="X1", "P"="X2"))

af_acat$FDR <- p.adjust(af_acat$P, method = "fdr")
af_acat$FWER <- p.adjust(af_acat$P, method="bonferroni")                                       

write.table(af_acat, file="~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MAGMA/Multiple_testing_corr/gProfiler_sets_GSA/Cauchy_combined_results.txt",
            sep = "\t", row.names = F, quote = F)
