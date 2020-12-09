########################################

## LD hub Volcano plot

########################################

library(readxl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

Pneumonia <- read.csv("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/LDSC_h2_rg/LDhub_pneumonia.csv", header = T)

ggplot(data = Pneumonia, 
       aes(x=rg, y=-log10(p), 
           colour = Category, label = trait2)) +
  geom_point( alpha = 0.8, size = 1.75) +
  labs(x = expression("Genetic correlation coefficient"), y = expression(paste("-log"[10], "P-value"))) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text_repel(aes(label=ifelse(p < 0.00027 & rg < 0, as.character(trait2), '')), force = 3.5, nudge_y = 1.1, nudge_x = -1, size = 2.5, segment.size = 0.2) +
  geom_text_repel(aes(label=ifelse(p < 0.00027 & rg > 0 & label=c("Insomina", "Obesity class 1", "Depressive Symptoms", "Neuroticism"), as.character(trait2), '')), force = 4, nudge_y = 1.2, nudge_x = 1.2, size = 2.5, segment.size = 0.2) +
  theme(axis.title = element_text(face="bold", size=12)) +
  geom_hline(yintercept = 3.69, linetype ="longdash") +
  ylim(c(0, 20)) +
  xlim(c(-0.6, 0.6))
