#############################

## Automated LCV of UKBB traits for pneumonia 2022

## William Reay (2022)

##############################


setwd("/home/control/data/users/william/2022_FinnGen_R6_meta/Meta_analysis_output/LCV")

##Load dependencies
suppressMessages(library(optparse))

##Specify command line inputs
option_list = list(
  make_option("--exposure_name", action="store", default=NA, type='character',
              help="The name of the exposure trait to be analysed [required]"),
  make_option("--outcome_name", action="store", default=NA, type='character',
              help="The name of the blood pressure trait to be analysed [required]"),
  make_option("--exposure_sumstats_file", action="store", default=NA, type='character',
              help="File name for exposure summary stats for trait to be analysed, gzip compression required"),
  make_option("--outcome_sumstats_file", action="store", default=NA, type='character',
              help="File name for outcome summary stats for trait to be analysed, gzip compression required"),
  make_option("--LD_score", action="store", default=NA, type='character',
              help="file with LD scores")
)
opt = parse_args(OptionParser(option_list=option_list))

#Load in munged summary statistics for FVC and trait to be analysed

cat("#########################")
cat("\n")
cat("Loading summary statistics for", paste(opt$outcome_name), "and", paste(opt$exposure_name))
cat("\n")
cat("#########################")
cat("\n")

#Load trait 1 data 
opt$exposure_name_df <- read.table(file = paste("../LDSR_UKBB_Neale_general/LDSR_Neale_UKBB/", opt$exposure_sumstats_file, sep=""),
                                   header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))

opt$exposure_name_df <- na.omit(opt$exposure_name_df)

#Load trait 2 data 
opt$outcome_name_df <- read.table(file = paste("../HapMap3_munged/", opt$outcome_sumstats_file, sep=""),
                                  header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))
opt$outcome_name_df <- na.omit(opt$outcome_name_df)

cat("\n")
cat("#########################")
cat("\n")
cat("GWAS data loaded")
cat("\n")
cat("#########################")
cat("\n")

cat("\n")
cat("#########################")
cat("\n")
cat("Loading LD scores")
cat("\n")
cat("#########################")
cat("\n")

#Load LD scores
LD_scores <- read.table(file = paste("../", opt$LD_score, sep=""),
                        header=TRUE,sep="\t",stringsAsFactors = FALSE, na.strings=c("", "NA"))


cat("\n")
cat("#########################")
cat("\n")
cat("Merging data frames")
cat("\n")
cat("#########################")
cat("\n")


#Merge such that SNPs are annotated with LD scores
Merged_df <- merge(LD_scores,opt$exposure_name_df,by="SNP")
Annotated_df_combined <- merge(Merged_df,opt$outcome_name_df,by="SNP")

#Sort by position 
Sorted_df <- Annotated_df_combined[order(Annotated_df_combined[,"CHR"],Annotated_df_combined[,"BP"]),]

Sorted_df$L2 <- as.numeric(Sorted_df$L2)

#Check if any mismatches
mismatch = which(Sorted_df$A1.x!=Sorted_df$A1.y,arr.ind=TRUE)
Sorted_df[mismatch,]$Z.y = Sorted_df[mismatch,]$Z.y*-1
Sorted_df[mismatch,]$A1.y = Sorted_df[mismatch,]$A1.x
Sorted_df[mismatch,]$A2.y = Sorted_df[mismatch,]$A2.x


cat("\n")
cat("#########################")
cat("\n")
cat("Constructing LCV model, output will be sent to text file",  paste(opt$exposure_name, "_LCV_", paste(opt$outcome_name), ".txt", sep = ""))
cat("\n")
cat("#########################")
cat("\n")


#Run LCV-need to setwd to directory containing LCV package
source("RunLCV.R")

file = paste('', opt$exposure_name, '_LCV_', paste(opt$outcome_name), '.txt', sep = "")

sink(file, append = TRUE)

LCV <- RunLCV(Sorted_df$L2,Sorted_df$Z.x, Sorted_df$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), pvalue=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, LCV$pval.gcpzero.2tailed, LCV$rho.est, LCV$rho.err)

#Sink output to text file

LCV

sink()

#Clear environment after run
rm(list = ls())