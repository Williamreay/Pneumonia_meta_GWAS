##############################

## LCV volcano plot

## William Reay (2020)

################################

library(readxl)
library(ggplot2)

LCV_figure <- read_excel("Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/LDSC_h2_rg/LCV_figure.xlsx")

ggplot(data = LCV_figure, aes(x=GCP, y=Z, colour=rg)) +
  ylim(0, 5) + xlim(-0.8, 0.8) +
  geom_hline(yintercept = 1.9666, linetype ="longdash") +
  geom_vline(xintercept=0.6, linetype="longdash") +
  geom_vline(xintercept=-0.6, linetype="longdash") +
  labs(x="Genetic causality proportion (GCP)", y="|Z score|", size="|Z|", colour="Genetic correlation") +
  theme_bw() +
  geom_point(aes(size=abs(Z))) +
  scale_colour_gradient2(low = "red",
                         mid = "white",
                         high = "blue") +
  geom_text(aes(label=ifelse(Z > 1.95 & Z < 2, as.character(Trait), '')), hjust=1.1, vjust=1.2, colour="black", size=3) +
  geom_text(aes(label=ifelse(Z > 3, as.character(Trait), '')), hjust=1.4, vjust=1.2, colour="black", size=3)
