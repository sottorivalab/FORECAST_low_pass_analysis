# Figure for recreating the KM curve for Figure 1

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(survival)
library(survminer)

# 250 patients outcome data
outcome = read.csv("refs/Fig1 FORECAST outcomes v4.csv")

colnames(outcome)[7:8] = c("Recurrence", "Recurrence_Months")

# Make a fit object for total events
fit = survfit(Surv(Recurrence_Months, Recurrence) ~ 1, data = outcome)

# Plot general KM curve
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = FALSE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               size = 0.6,
               censor.size = 3,
               ylab = paste0("Fraction Recurrence Free"),
               xlab = "Time (Years)",
               xscale = 12,
               break.x.by = 12*4,
               tables.theme = theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    plot.title = element_blank()),
               palette = c("black"))
ggsave(plot = print(p), "results/9_outcome_analysis/Figure_1_KM_curve_250_patients.pdf", width = 5, height = 3)
