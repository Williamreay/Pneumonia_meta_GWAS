################################

## Effect of PES on antibody response to antigens amongst seropositive individuals

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
  make_option("--antigen_name", action="store", default=NA, type='character',
              help="Name of antigen to be correlated with genetic score amongst seropositive individuals  [required"),
  make_option("--antigen_name_col_id", action="store", default=NA, type='character',
              help="Column id for the antigen [required"),
  make_option("--score_name", action="store", default=NA, type='character',
              help="Name of PGS or PES [required"),
  make_option("--score_col_id", action="store", default=NA, type='character',
              help="Column id for the PGS or PES [required")
)
opt = parse_args(OptionParser(option_list=option_list))

print(opt)

cat("\n")
cat("#########################")
cat("\n")
cat("Constructing model amongst",  paste(opt$infection_name, ' seropositive individuals that correlates ', opt$score_name, ' with ', opt$antigen_name, sep = ""))
cat("\n")
cat("#########################")
cat("\n")

## Read in RDS

Merged_antigen_data <- readRDS(file = "Cleaned_antigen_data.rds")


# MFI_PES_PRS <- function(Seropositivity_col, Antigen, df, score) {
#   var1 <- enquo(Seropositivity_col)
#   seropos <- df %>% filter(!! var1 == 1)
#   seropos_df <- as.data.frame(seropos)
#   seropos_df$Scaled_antigen <- as.numeric(scale(seropos_df$Antigen))
#   seropos_df$score <- as.numeric(scale(seropos_df$score))
#   ##Test association between score and antigen MFI amongst seropositive individuals
#   mod <- lm(Scaled_antigen ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
#             PC8 + PC9 + PC10 + Batch + Freeze_thaw + Spill_over + score, data = seropos_df)
#   return(mod)
# }

# seropos_col <- "f.23050.0.0"
# antigen_col <- "f.23000.0.0"
# score_col <- "Bile_acid_PES"
# df <- Merged_antigen_data

                                              
MFI_PES_PRS <- function(seropos_col, antigen_col, score_col, df) {
  seropos_df <- df[df[[seropos_col]] == 1,]
  seropos_df$scaled_antigen <- as.numeric(scale(seropos_df[[antigen_col]]))
  seropos_df$score <- as.numeric(scale(seropos_df[[score_col]]))
  mod <- lm(scaled_antigen ~ Sex + Age + Age2 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
              PC8 + PC9 + PC10 + Batch + Freeze_thaw + Spill_over + score, data = seropos_df)
  return(summary(mod))
}


Run_lm <- MFI_PES_PRS(opt$infection_name_col_id, opt$antigen_name_col_id, opt$score_col_id, Merged_antigen_data)

Output <- as.data.frame(Run_lm$coefficients)[18, 1:4]

write.table(Output, file = paste("Antigen_PES_PRS_output/", opt$infection_name, "_", 
                                 opt$antigen_name, "_", opt$score_name, ".txt", sep=""),
            sep = "\t", row.names = F, quote = F)

#Clear environment after run
rm(list = ls())
