#####################################

## LAVA local genetic correlation at genome-wide sig pneumonia loci

## V2G or closest gene

## MHC treated as single region

#####################################

library(LAVA)

setwd("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/LAVA_local_rg_psych_pheno/")

input = process.input(input.info.file = "Local_rg_GW_sig_input.txt",
                       ref.prefix = "./g1000_eur",
                       phenos = c("Pneumonia", "ADHD", "MDD", "PTSD"),
                       sample.overlap.file = NULL)

loci = read.loci("GW_sig_locus.txt")

## Iteratively estimate univariable local h2 for each locus

## MHC

locus =  process.locus(loci[1,], input)

MHC_univar_h2 <- run.univ(locus)

write.table(MHC_univar_h2, file="MHC_univariable_h2.txt", sep = "\t", row.names = F, quote = F)

MHC_bivar <- run.bivar(locus)

write.table(MHC_bivar, file="MHC_bivariate_corr_pneumonia_psych.txt", sep = "\t", row.names = F, quote = F)

## Chr 11 (mucin dense) locus

locus =  process.locus(loci[2,], input)

Chr_11_univar_h2 <- run.univ(locus)

write.table(Chr_11_univar_h2, file="Univar_h2_chr_11.txt", sep = "\t", row.names = F, quote = F)

Chr_11_bivar <- run.bivar(locus)

write.table(Chr_11_bivar, file="Bivar_rg_chr_11.txt", sep = "\t", row.names = F, quote = F)

## Chr 12 locus

locus =  process.locus(loci[3,], input)

Chr_12_univar_h2 <- run.univ(locus)

write.table(Chr_12_univar_h2, file="Univar_h2_chr_12.txt", sep = "\t", row.names = F, quote = F)

## Chr 5 locus

locus =  process.locus(loci[4,], input)

Chr_5_univar_h2 <- run.univ(locus)

write.table(Chr_5_univar_h2, file="Univar_h2_chr_5.txt", sep = "\t", row.names = F, quote = F)

## Chr 1 locus

locus =  process.locus(loci[5,], input)

Chr_1_univar_h2 <- run.univ(locus)

write.table(Chr_1_univar_h2, file="Univar_h2_chr_1.txt", sep = "\t", row.names = F, quote = F)


