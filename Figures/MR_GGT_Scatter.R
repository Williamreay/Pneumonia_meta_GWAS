## Custom scatter plot for GGT on pneumonia susceptibility


## Get Egger Slope

Slope <- mr_egger_regression(GGT_pneu_harm$beta.exposure, GGT_pneu_harm$beta.outcome, GGT_pneu_harm$se.exposure, GGT_pneu_harm$se.outcome, default_parameters())

Univariable_GGT$intercept <- c(0, 0, 0, 0, -0.0008042)

Univariable_GGT$METHOD <- as.factor(c("IVW-MRE", "IVW-FE", "Weighted Median", "Weighted Mode", "MR-Egger"))

ggplot(data=GGT_pneu_harm, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(aes(text=paste("SNP:", SNP))) +
  theme_bw() +
  geom_abline(data=Univariable_GGT, aes(intercept=intercept, slope=b, colour=METHOD), show.legend=TRUE) +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
  labs(colour="MR Test", x="SNP effect on GGT", y="SNP effect on pneumonia") +
  theme(legend.position="right", legend.direction="vertical") +
  guides(colour=ggplot2::guide_legend(ncol=1)) +
  xlim(0, 0.11)

## Triglycerides

ggplot(data=TG_pneumonia_harmonised, aes(x=beta.exposure, y=beta.outcome)) +
  geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
  geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
  geom_point(aes(text=paste("SNP:", SNP))) +
  theme_bw() +
  geom_abline(data=TG_to_pneumonia, aes(intercept=intercept, slope=b, colour=METHODS), show.legend=TRUE) +
  scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
  labs(colour="MR Test", x="SNP effect on triglycerides", y="SNP effect on pneumonia") +
  theme(legend.position="right", legend.direction="vertical") +
  guides(colour=ggplot2::guide_legend(ncol=1)) +
  xlim(0, 0.15) + ylim(-0.06, 0.06)

## Get Egger Slope

Slope <- mr_egger_regression(TG_pneumonia_harmonised$beta.exposure, TG_pneumonia_harmonised$beta.outcome, TG_pneumonia_harmonised$se.exposure, TG_pneumonia_harmonised$se.outcome, default_parameters())

TG_to_pneumonia$intercept <- c(0, 0, 0, 0, 0.0009949224)
