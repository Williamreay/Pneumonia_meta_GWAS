##############################

## Cauchy aggregation of rare variants at gene level

## William Reay (2020)

###############################

library(dplyr)

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/Rare_var/Cauchy/")

## Load in all rare variants and their P value

Rare_var_p_val <- read.table("Filtered_rare_var.bed", header = T)

## Load in all annotated rare variants and drop intergenic variants

All_rare_variants_RAW <- read.table("Unfiltered_rare_variants//Pneumonia_gene_anno_rare_var.variant_function", header = T)

All_rare_variants_final <- All_rare_variants_RAW %>% filter(Annotation != 'intergenic')

## Merge rare variants with p values

Merged_all_rare_variants <- merge(All_rare_variants_final, Rare_var_p_val, by = "SNP" )

## Load NCBI 37 genes

Gene_anno <- read.table("NCBI37_hg19_genes.bed", header = T)

## Merge with NCBI protein coding genes

Merged_all_rare_var_gene_anno <- merge(Merged_all_rare_variants, Gene_anno, by ="Gene")

genes <- sort(unique(Merged_all_rare_var_gene_anno$Gene))
gene_ps <- list()

for (i in genes) { gene_ps[[i]] <- c(Merged_all_rare_var_gene_anno[Merged_all_rare_var_gene_anno$Gene==i,][12])}

genes_acat <- list()
for (i in genes) { genes_acat[[i]] <- ACAT(unlist(gene_ps[[i]], use.names=FALSE))}

af_acat <- data.frame(cbind(as.vector(genes), as.vector(unlist(genes_acat))))

af_acat <- rename(af_acat, c("Gene"="X1", "P"="X2"))

af_acat$FDR <- p.adjust(af_acat$P, method = "fdr")
af_acat$FWER <- p.adjust(af_acat$P, method="bonferroni")

## Write ACAT results

write.table(af_acat, file="Output/All_rare_var_ACAT_no_intergenic.txt",
          row.names = F, quote = F, sep = "\t")

## Write Merged output

write.table(Merged_all_rare_var_gene_anno, file="Output/Rare_variant_annotation.txt",
            sep = "\t", row.names = F, quote = F)
#'
#' Aggregated Cauchy Assocaition Test
#'
#' A p-value combination method using the Cauchy distribution.
#'
#'
#'
#' @param Weights a numeric vector of non-negative weights for the combined p-values. When it is NULL, the equal weights are used.
#' @param Pvals a numeric vector of p-values to be combined by ACAT.
#' @return p-value of ACAT.
#' @author Yaowu Liu
#' @examples p.values<-c(2e-02,4e-04,0.2,0.1,0.8)
#' @examples ACAT(Pvals=p.values)
#' @export
ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }
  
  
  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}
