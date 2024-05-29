# Run 11_outcome_analysis.R until line 106

# Start with Time to Recurrence
out_var = "Recurrence"
out_tim = paste0(out_var,"_Months")

# Data transformation and formatting
collected$Lossness                = scaleQuant(log(collected$Lossness))
collected$Spearman                = scaleQuant(collected$Spearman)
collected$High_events             = collected$Total_Events > median(collected$Total_Events)
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 5:2)

# List covariates
covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", "Lossness", "High_events", "Spearman")

# Now, using the same covariates, perform a multivariate test
res.cox = coxph(as.formula(paste0("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected)

# Collect the results for a table 
sum.res.cox = summary(res.cox)

# Format
res.table = signif(cbind(sum.res.cox[["conf.int"]][,c(1,3,4)], sum.res.cox[["coefficients"]][,"Pr(>|z|)"]), digits = 3)

# Name rows and columns
colnames(res.table) = c("HR", "Lower", "Upper", "p")
rownames(res.table) = c("Samples", "PSA_20+", "Gleason_4", "Gleason_3", "Gleason_2", "T3+", "N1+", "Lossness", "High_Events", "Spearman")

# Format R object
res.table = as.data.frame(res.table)

# Make more human readable
res.table$Reference = c("0", "PSA_less_than_20", 
                        "Gleason_5", "Gleason_5", "Gleason_5", 
                        "Not_T3", "N0", "0", "Low_Events", "0")

# Output result 
write.csv(res.table, file = "results/9_outcome_analysis/All_the_best_CPH_recurrence.csv", quote = F)

# Change this back for plotting ordering
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 2:5)

# Use survminer to plot the results
p = ggforest(res.cox)

# Save ggforest
ggsave(plot = p, "results/9_outcome_analysis/All_the_best_CPH_recurrence.pdf", width = 8, height = 7)

#########################################################################################################################################

# Read in the sidedness results
phylosig_side = read.csv("results/8_collect_metrics/FORECAST_phylosig_sidedness_analysis.csv", 
                         stringsAsFactors = F)

# Filter for successful assessment
phylosig_side = phylosig_side[!is.na(phylosig_side$p),]

# Add sidedness lambda metrics
collected = merge(collected, phylosig_side[,c("Patient", "lambda")], by = "Patient", all.x = T)

# Add category
collected$phylosig = collected$lambda > 0.8

# Test phylosig analysis in a multivariate
covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", "phylosig", "Spearman")

# Remove NAs and create an alt object
collected_sub = collected[!is.na(collected$phylosig),]

# Make 5 reference
collected_sub$Gleason_grade_group     = factor(collected_sub$Gleason_grade_group, levels = 5:2)

# Now perform a multivariate test
res.cox = coxph(as.formula(paste("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected_sub)

# For plotting
collected_sub$Gleason_grade_group     = factor(collected_sub$Gleason_grade_group, levels = 2:5)

# Use survminer to plot the results
p = ggforest(res.cox)

# Save ggforest
ggsave(plot = p, "results/9_outcome_analysis/Phylosig_CPH_recurrence.pdf", width = 8, height = 7)

#########################################################################################################################################

# Now we do Time to Metastasis
out_var = "Mets"
out_tim = paste0(out_var,"_Months")

# Data transformation and formatting
collected$High_spearman           = collected$Spearman > quantile(collected$Spearman, probs = 2/3)
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 5:2)

# What are the covariates?
covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", 
               "Lossness", 
               "High_events",
               "High_spearman", 
               "MYC_FGFR1")

# Now, using the same covariates, perform a multivariate test
res.cox = coxph(as.formula(paste0("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected)

# Reorder for plotting
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 2:5)

# Use survminer to plot the results
p = ggforest(res.cox)

# Output ggforest
ggsave(plot = p, "results/9_outcome_analysis/All_the_best_CPH_mets.pdf", width = 8, height = 7)
