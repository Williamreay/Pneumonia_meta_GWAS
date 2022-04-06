## Custom scatter plot for CRP on pneumonia susceptibility


## Get Egger Slope

Slope <- mr_egger_regression(CRP_pneu_harm$beta.exposure, CRP_pneu_harm$beta.outcome, CRP_pneu_harm$se.exposure, CRP_pneu_harm$se.outcome, default_parameters())

Univariable_CRP$intercept <- c(0, 0, 0, 0, 0.002257)

Univariable_CRP$METHOD <- as.factor(c("IVW-MRE", "IVW-FE", "Weighted Median", "Weighted Mode", "MR-Egger"))

ggplot(data=CRP_pneu_harm, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(aes(text=paste("SNP:", SNP))) +
  theme_bw() +
  geom_abline(data=Univariable_CRP, aes(intercept=intercept, slope=b, colour=METHOD), show.legend=TRUE) +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
  labs(colour="MR Test", x="SNP effect on CRP", y="SNP effect on pneumonia") +
  theme(legend.position="right", legend.direction="vertical") +
  guides(colour=ggplot2::guide_legend(ncol=1)) +
  xlim(0, 0.11)

