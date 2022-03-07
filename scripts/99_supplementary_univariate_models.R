# Where's the target data?
target_repo = "target_repo/"

##########################################################
# Perform recurrence part of 99_outcome_analysis.R first #
##########################################################

collected$High_PGA = collected$mPGA > median(collected$mPGA)

# Make a fit object for total events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_PGA")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
           pval = TRUE, conf.int = FALSE,
           legend = "none",
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = c(1,1), # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           tables.y.text = FALSE,
           title = "Mean PGA",
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
ggsave(plot = print(p), "results/11_outcome_analysis/PGA_median_recurrence.png", width = 4, height = 4)

##########################################################

collected$High_PGA = collected$mPGA >= 7.49/100

# Make a fit object for total events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_PGA")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "Mean PGA (>= 7.49%)",
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
ggsave(plot = print(p), "results/11_outcome_analysis/PGA_Lalonde_split_recurrence.png", width = 4, height = 4)

##########################################################

collected$High_max_PGA = collected$max_PGA > median(collected$max_PGA)

# Make a fit object for total events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_max_PGA")), data = collected)

# Plot the results of the total event threshold
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "Max PGA",
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
ggsave(plot = print(p), "results/11_outcome_analysis/Max_PGA_median_recurrence.png", width = 4, height = 4)

##########################################################

fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYC")), data = collected)

p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "MYC amplification",
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
ggsave(plot = print(p), "results/11_outcome_analysis/MYC_amp_recurrence.png", width = 4, height = 4)

##########################################################

fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ ATM")), data = collected)

p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "ATM mutation",
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
ggsave(plot = print(p), "results/11_outcome_analysis/ATM_mutation_recurrence.png", width = 4, height = 4)

##########################################################

fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ Subclonal_Mut_Status")), data = collected)

p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "Subclonal Driver Mutation",
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
ggsave(plot = print(p), "results/11_outcome_analysis/Subclonal_driver_mutation_recurrence.png", width = 4, height = 4)

##########################################################

# Read in the heatmap 
driver = read.csv(paste0(target_repo,"results/mutation_calling/Per_patient_mutation_status_heatmap_data.csv"), 
                  stringsAsFactors = F, row.names = "X")

driver_number = apply(driver, 2, function(i) length(i[na.omit(i)!=0]))

driver_number = data.frame(Patient = names(driver_number), Number_Panel_Mutations = driver_number)

rownames(driver_number) = NULL

collected = merge(collected, driver_number, by = "Patient")

collected$Number_Panel_Mutations[is.na(collected$ATM)] = NA

# The best separation
collected$High_panel_number = collected$Number_Panel_Mutations > 1

fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_panel_number")), data = collected)

p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               title = "Number of panel mutations (>1)",
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
ggsave(plot = print(p), "results/11_outcome_analysis/Number_panel_mutations_recurrence.png", width = 4, height = 4)
