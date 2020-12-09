###########################
#Drug target-gene interaction enrichment for candidate PES pathways

#William Reay - November 2019 (william.reay@uon.edu.au)
###########################

#Load dependencies

list.of.packages <- c("optparse", "WebGestaltR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(optparse))
library(WebGestaltR)

#Specify command line inputs
option_list = list(
  make_option("--geneset_id", action="store", default=NA, type='character',
              help="The name of the geneset of interest, case sensitive [required]"),
  make_option("--phenotype", action="store", default=NA, type='character',
              help="Name of phenotype of interest"),
  make_option("--path",action="store", default=NA, type='character',
              help="File path to list of genes")
)
opt = parse_args(OptionParser(option_list=option_list))

## Set working directory

setwd(opt$path)

## Perform overrepresentation analysis, FDR correction

ORA_DrugBank <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens", 
                            enrichDatabase = "drug_DrugBank", 
                            enrichDatabaseType = "genesymbol", 
                            interestGeneFile = paste("", opt$geneset_id, "/", opt$geneset_id, "_genes.txt", sep= ""), 
                            interestGeneType = "genesymbol", 
                            referenceSet = "genome", minNum = 3, 
                            sigMethod = "fdr", fdrMethod = "BH", fdrThr = 0.05,
                            outputDirectory = paste("",opt$path, "", opt$geneset_id,"/", sep=""), projectName = paste("DrugBank_ORA_",opt$geneset_id, "_", opt$phenotype, sep=""))

rm(list = ls())