#################################

## SNP h2 plots

## William Reay (2020)

#################################

library(ggplot2)
library(readxl)

LDSC_h2_liability <- read_excel("~/Desktop/23andMe_pneumonia/FINAL_IVW_meta_analysis/LDSC_h2_rg/LDSC_h2_liability.xlsx")

h2_liab <- ggplot(LDSC_h2_liability, aes(x=Source, y=h2_liab, fill=Source)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=h2_liab-h2_liab_se, ymax=h2_liab+h2_liab_se), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  labs(y=expression(paste(h^2, " liability scale (+/- SE)"))) +
  theme( axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
