library(data.table)
library(ggplot2)
library(dplyr)

IN <- fread("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/MR_pheWAS_drug_targets/TNFRSF1A_ieugwas_MR_PheWAS.txt")

IN <- IN %>% filter(MR_FWER < 0.05)

IN$Category <- c("Haemtological", "Haemtological", "Haemtological",
                 "Haemtological", "Haemtological", "Haemtological",
                 "Haemtological", "Inflammatory disease",
                 "Haemtological", "Haemtological", "Haemtological", 
                 "Haemtological", "Inflammatory disease", 
                 "Inflammatory disease", "Haemtological", "Haemtological",
                 "Haemtological", "Haemtological", "Inflammatory disease")

IN$U_CI <- IN$MR_beta + (1.96*IN$MR_SE)
IN$L_CI <- IN$MR_beta - (1.96*IN$MR_SE)

IN$MR_Z <- (IN$MR_beta/IN$MR_SE)

IN$Flipped <- IN$MR_Z*-1

IN <- IN %>% filter(id != "ieu-a-1025")

IN <- IN %>% filter(trait != "Monocyte cell count" & trait != "Lymphocyte cell count" &
                      trait != "Eosinophil cell count")

ggplot(data = IN, aes(trait, Flipped, fill = Category, colour = Category)) +
  geom_bar(stat = "identity") + coord_flip() +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("Effect of genetically proxied TNFRSF1A inhibition (Z score)") +
  scale_fill_manual(values = c("steelblue1", "grey")) +
  scale_colour_manual(values= c("black", "black")) +
  xlab(" ")

