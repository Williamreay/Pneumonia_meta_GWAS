################################

## Effect of PES on seropositivity

## William Reay (2020)

#################################

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

setwd("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/UKBB_pneumonia/PRS/Antigen_data/")


##Specify command line inputs
option_list = list(
  make_option("--infection_name", action="store", default=NA, type='character',
              help="The name of infection that will be analysed to identify seropositive UKBB partcipants [required]"),
  make_option("--infection_name_col_id", action="store", default=NA, type='character',
              help="Column id for the infection [required"),
  make_option("--score_name", action="store", default=NA, type='character',
              help="Name of PGS or PES [required"),
  make_option("--score_col_id", action="store", default=NA, type='character',
              help="Column id for the PGS or PES [required")
)
opt = parse_args(OptionParser(option_list=option_list))


cat("\n")
cat("#########################")
cat("\n")
cat("Constructing to test association between",  paste(opt$score_name, ' and seropositivity for ', opt$infection_name, sep = ""))
cat("\n")
cat("#########################")
cat("\n")

## Read in RDS

Merged_antigen_data <- readRDS(file = "Cleaned_antigen_data.rds")

Seropos_PES_PRS <- function(seropos_col, score_col, df) {
  df$score <- as.numeric(scale(df[[score_col]]))
  df$seropos <- df[[seropos_col]]
  mod <- glm(seropos ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + Batch + Freeze_thaw + Spill_over + score, family = "binomial", data = df)
  return(summary(mod))
}

Run_glm <- Seropos_PES_PRS(opt$infection_name_col_id, opt$score_col_id, Merged_antigen_data)

Output <- as.data.frame(Run_glm$coefficients)[18, 1:4]

write.table(Output, file = paste("Seropos_output/", opt$infection_name, "_", 
                                 opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = F, quote = F)

#Clear environment after run
rm(list = ls())
