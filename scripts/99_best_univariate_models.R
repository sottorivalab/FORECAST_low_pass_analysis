##########################################################
# Perform recurrence part of 11_outcome_analysis.R first #
##########################################################

# Apply thresholds
collected$High_Tot_Events = collected$Total_Events > median(collected$Total_Events)
collected$High_spearman   = collected$Spearman > quantile(collected$Spearman, probs = 2/3)

# Make a fit object for total events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_Tot_Events")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           legend = "none",
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = c(1,1), # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           tables.y.text = FALSE,
           title = "Tree Events",
           ylab = paste0("Fraction ", out_var," Free"),
           xlab = "Time (Years)",
           xscale = 12,
           break.x.by = 12*4,
           tables.theme = theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                plot.title = element_blank()),
           ggtheme = theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)), # Change ggplot2 theme
           palette = c("#929591", "#e50000"))
ggsave(plot = print(p), "results/9_outcome_analysis/Events_median_recurrence.png", width = 4, height = 4)

# Make a fit object for the Spearman split
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_spearman")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "Spearman's rho",
               ylab = paste0("Fraction ", out_var," Free"),
               xlab = "Time (Years)",
               xscale = 12,
               break.x.by = 12*4,
               tables.theme = theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    plot.title = element_blank()),
               ggtheme = theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)), # Change ggplot2 theme
               palette = c("#929591", "#e50000"))
ggsave(plot = print(p), "results/9_outcome_analysis/Spearman_tertile_recurrence.png", width = 4, height = 4)

###############################
# Phylosig sidedness analysis #
###############################

# Read in the sidedness results
phylosig_side = read.csv("results/8_collect_metrics/FORECAST_phylosig_sidedness_analysis.csv", 
                         stringsAsFactors = F)

# Filter for successful assessment
phylosig_side = phylosig_side[!is.na(phylosig_side$p),]

# Add sidedness lambda metrics
collected = merge(collected, phylosig_side[,c("Patient", "lambda")], by = "Patient", all.x = T)

# Add category
collected$phylosig = collected$lambda > 0.8

# Make a fit object
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ phylosig")), data = collected)

# Plot the results of phylosig sidedness
p = ggsurvplot( fit,
                pval = TRUE, conf.int = FALSE,
                legend = "none",
                risk.table = TRUE, # Add risk table
                risk.table.col = "strata", # Change risk table color by groups
                linetype = c(1,1), # Change line type by groups
                surv.median.line = "hv", # Specify median survival
                tables.y.text = FALSE,
                title = "Phylosig sidedness",
                ylab = y_label,
                xlab = "Time (Years)",
                xscale = 12,
                break.x.by = 12*4,
                tables.theme = theme(axis.title.x=element_blank(),
                                     axis.text.x=element_blank(),
                                     axis.ticks.x=element_blank(),
                                     plot.title = element_blank()),
                ggtheme = theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)), # Change ggplot2 theme
                palette = c("#929591", "#e50000"))
ggsave(plot = print(p), paste0("results/9_outcome_analysis/Phylosig_sidedness_recurrence.png"), 
       width = 4, height = 4)

##########################################################
# Perform metastasis part of 11_outcome_analysis.R first #
##########################################################

# Make a fit object for the % of subclonal events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYC_FGFR1")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "MYC/FGFR1 amplification",
               ylab = paste0("Fraction ", ifelse(out_var=="Mets", "Metastasis", out_var)," Free"),
               xlab = "Time (Years)",
               xscale = 12,
               break.x.by = 12*4,
               tables.theme = theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    plot.title = element_blank()),
               ggtheme = theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)), # Change ggplot2 theme
               palette = c("#929591", "#e50000"))
ggsave(plot = print(p), "results/9_outcome_analysis/Myc_Fgfr1_metastasis.png", width = 4, height = 4)
