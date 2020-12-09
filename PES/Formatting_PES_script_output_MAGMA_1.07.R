####################

## Load .sets.out or .gsa.out files into R and perform multiple testing correction

## Works for v1.07 filenames

## William Reay (2020)

## Check if dependencies are present, if not, install them

list.of.packages <- c("optparse", "data.table", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Load dependencies
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))

## Specify command line inputs
option_list = list(
  make_option("--phenotype_name", action="store", default=NA, type='character',
              help="The name of the trait to be analysed [required]"),
  make_option("--path", action="store", default = NA, type='character',
              help="Path to PES script output"),
  make_option("--geneset_file_prefix",action="store", default=NA, type='character',
              help="Prefix of .sets.out file, no . needed")
)
opt = parse_args(OptionParser(option_list=option_list))

setwd(paste(opt$path))

## Import output file for all SNPs

All_SNPs <- fread(file=paste("",opt$geneset_file_prefix,"allSNPs.gsa.out", sep=""), skip=3,
                  stringsAsFactors = F, data.table = F)

## Make a column which denotes the P value threshold

All_SNPs$P_threshold <- 1

## Import output file for P < 0.5

SNPs_below_0.5 <- fread(file=paste("",opt$geneset_file_prefix,"_0.5_SNPs.gsa.out", sep=""), skip=3,
                        stringsAsFactors = F, data.table = F)

## Make a column which denotes the P value threshold

SNPs_below_0.5$P_threshold <- 0.5

## Import output file for P < 0.05

SNPs_below_0.05 <- fread(file=paste("",opt$geneset_file_prefix,"_0.05_SNPs.gsa.out", sep =""), skip=3,
                         stringsAsFactors = F, data.table = F)

## Make a column which denotes the P value threshold

SNPs_below_0.05$P_threshold <- 0.05

## Import output file for P < 0.005

SNPs_below_0.005 <- fread(file=paste("",opt$geneset_file_prefix,"_0.005_SNPs.gsa.out", sep = ""), skip=3,
                          stringsAsFactors = F, data.table = F)

## Make a column which denotes the P value threshold

SNPs_below_0.005$P_threshold <- 0.005

## Vertically concatenate files 

Combined_gsa <- bind_rows(All_SNPs, SNPs_below_0.5, SNPs_below_0.05, SNPs_below_0.005)

## Select sets with > 5 genes

Combined_gsa <- Combined_gsa %>% filter(Combined_gsa$NGENES > 5)

## Apply multiple testing correction (FDR)

Combined_gsa$FDR <- p.adjust(Combined_gsa$P, method="fdr")

## Apply Bonferroni correction

Combined_gsa$FWER_P <- p.adjust(Combined_gsa$P, method="bonferroni")

## Write column for sets which survive multiple testing correction FDR < 0.05

Combined_gsa$Multiple_testing <- ifelse(Combined_gsa$FDR < 0.05, "YES", "NO")

Combined_gsa$Multiple_testing_FWER <- ifelse(Combined_gsa$FWER_P < 0.05, "YES", "NO")

cat("#########################")
cat("\n")
cat("Printing output for sets with > 5 overlapping genes")
cat("\n")
cat("#########################")
cat("\n")

## Write output to a .csv

write.csv(Combined_gsa, file=paste(opt$path,"",opt$phenotype_name, "_all_Tclin_sets.csv", sep=""), quote = F, row.names = F)

rm(list = ls())
