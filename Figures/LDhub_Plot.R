###################################

## LDSC LD hub plot

## William Reay (2020)

####################################

library(ggplot2)
library(readxl)


LDSC_FP_input <- read_excel("Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/LDSC_h2_rg/LDSC_FP_input.xlsx")

LDSC_FP_input$Urg <- LDSC_FP_input$rg + LDSC_FP_input$se
LDSC_FP_input$Lrg <- LDSC_FP_input$rg - LDSC_FP_input$se

FP <- ggplot(data = LDSC_FP_input, aes(x=trait2, y=rg, ymin=Urg, ymax=Lrg, colour=Category)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +
  facet_wrap(~Category,strip.position="top",nrow=4,scales = "free_y") +
  coord_flip() +
  xlab("") + ylab("Genetic correlation estimate") +
  theme_bw() +
  theme(legend.position = "none")
  labs(colour="Category")
