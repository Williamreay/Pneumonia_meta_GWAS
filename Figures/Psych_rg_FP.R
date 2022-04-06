## Psych rg forest plot

library(ggplot2)
library(readxl)
library(viridis)

DF <- read_excel("~/Desktop/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Common_var_results/LDSR/PGC_psych_LDSR_with_and_without_smoking_conditioning/PGC_rg_results.xlsx")

DF$UCI <- DF$rg + (1.96*DF$se)
DF$LCI <- DF$rg - (1.96*DF$se)

FP <- ggplot(data = DF, aes(x=Trait, y=rg, ymin=LCI, ymax=UCI, colour=Trait)) +
  geom_pointrange() +
  geom_hline(yintercept=0, lty=2) +
  coord_flip() +
  ylab("Bivariate genome-wide genetic correlation estimate") +
  theme_bw() +
  theme(legend.position = "null", axis.title.y = element_blank()) +
  facet_wrap(~GWAS,strip.position="top",nrow=2,scales = "free_y")

FP + scale_color_viridis(discrete = T)

