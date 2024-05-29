# Work looking at comparing Gleason Morisita to the Spearman metric

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(survival)
library(survminer)

# What outcome?
out_tim = "Recurrence_Months"
out_var = "Recurrence"

# Read in the genomic metrics
patient_mets = read.table("results/8_collect_metrics/FORECAST_genomic_metrics.txt", sep = "\t", header = T, stringsAsFactors = F)

# Get the gleason morisita
gleason_morisita = read.csv("refs/GleasonMorisita_Patients.csv", header = F, stringsAsFactors = F)

# Make the column names readable
colnames(gleason_morisita) = c("Patient", "Gleason_Morisita")

# Merge it together
patient_mets = merge(patient_mets, gleason_morisita)

# Take only those with enough lpWGS samples
patient_mets = patient_mets[patient_mets$Samples > 2,]

# Linear model of PGA and GM
lm = summary(lm(Gleason_Morisita ~ mPGA, patient_mets[patient_mets$Gleason_Morisita!=0,]))

# What's the p value?
p_value = signif(lm$coefficients[2,4], digits = 2)

# Save that plot
ggplot(patient_mets[patient_mets$Gleason_Morisita!=0,], aes(x = mPGA, y = Gleason_Morisita)) + 
   annotate(geom="text", x=0.4, y=0.2, label=paste0("p=",format(p_value, scientific = FALSE)), size = 5) +
   ylab("Gleason Morisita") +
   geom_point() + geom_smooth(method = "lm", col = "#0165fc", se = T)
ggsave("results/8_collect_metrics/Gleason_Morisita_mPGA.png", height = 3, width = 4)

# Linear model again for Spearman and GM this time
lm = summary(lm(Gleason_Morisita ~ Spearman, patient_mets[patient_mets$Gleason_Morisita!=0,]))

# P value?
p_value = signif(lm$coefficients[2,4], digits = 2)

# Make that plot
ggplot(patient_mets[patient_mets$Gleason_Morisita!=0,], aes(x = Spearman, y = Gleason_Morisita)) + 
   annotate(geom="text", x=0.4, y=0.5, label=paste0("p=",format(p_value, scientific = FALSE)), size = 5) +
   ylab("Gleason Morisita") +
   geom_point() + geom_smooth(method = "lm", col = "#0165fc", se = T)
ggsave("results/8_collect_metrics/Gleason_Morisita_Spearman.png", height = 3, width = 4)

# Read in the outcome data
outcome = read.csv("refs/FORECAST_outcome_data_v2.csv", stringsAsFactors = F)

# Read in the data characteristics
characteristics = read.csv("refs/FORECAST_patient_characteristics.csv", stringsAsFactors = F)

# Read in the gleason grade group
digipath_gleason = read.csv("refs/GleasonGradeGroupPatient.csv", 
                            header = F, stringsAsFactors = F)
colnames(digipath_gleason) = c("Patient", "Digital_Gleason")
digipath_gleason$Digital_Gleason = floor(digipath_gleason$Digital_Gleason)

# Read in the rescoring pathologist scores
gleason_rescore = read.csv("refs/GleasonRescore_patients.csv", 
                           header = T, stringsAsFactors = F)

# Make columns more readable
colnames(outcome) = c("Patient", "Recurrence", "Recurrence_Months", "Death", "Death_Months", "Mets",
                      "Mets_Months", "Wide_Mets", "Wide_Met_Months")

# Rename the characteristics column
colnames(characteristics)[1] = "Patient"

# Add IM to beginning
outcome$Patient         = ifelse(grepl("DNT", outcome$Patient), 
                                 yes = outcome$Patient, no = paste0("IM",outcome$Patient))
characteristics$Patient = ifelse(grepl("DNT", characteristics$Patient),
                                 yes = characteristics$Patient, no = paste0("IM",characteristics$Patient))

# Merge by patient id
collected = merge(outcome, characteristics, by = "Patient")
collected = merge(collected, patient_mets, by = "Patient")
collected = merge(collected, digipath_gleason, by = "Patient")
collected = merge(collected, gleason_rescore[,c("Patient", "Gleason_Grade_Rescore")], by = "Patient", all.x = T)

# Remove anything that doesn't have enough samples
collected = collected[collected$Samples > 2,]

# Double category
collected$Morpho_Genome_het = apply(collected, 1, function(r) paste0(r["Spearman"]>median(collected$Spearman),"_",r["Gleason_Morisita"]>median(collected$Gleason_Morisita)))

# More informative names? Kind of...
collected$Morpho_Genome_het[collected$Morpho_Genome_het=="FALSE_FALSE"] = c("G_Hom_M_Hom")
collected$Morpho_Genome_het[collected$Morpho_Genome_het=="FALSE_TRUE"] = c("G_Hom_M_Het")
collected$Morpho_Genome_het[collected$Morpho_Genome_het=="TRUE_TRUE"] = c("G_Het_M_Het")
collected$Morpho_Genome_het[collected$Morpho_Genome_het=="TRUE_FALSE"] = c("G_Het_M_Hom")

# Calculate fit for quartiles
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,")~ Morpho_Genome_het")), data = collected)

# Plot quartiles of spearman/gm heterogeneity
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               risk.table = FALSE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1,1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               legend = "none",
               ylab = paste0("Fraction ",out_var," Free"),
               xlab = "Time (Months)",
               palette = c("#e50000", "#929591", "#6b8ba4", "#363737"))
ggsave(plot = print(p), "results/9_outcome_analysis/Spearman_GM_quartiles_survival.png", width = 4, height = 3)

# Is it doubly heterogeneous?
collected$Morpho_Genome_bad = collected$Morpho_Genome_het == "G_Het_M_Het"

# Calculate fit for most het vs. rest
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,")~ Morpho_Genome_bad")), data = collected)

# Plot double bad result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           palette = c("#0343df", "#e50000"),
           ggtheme = theme_bw())

# Do another univariate but this time with lower tertile
collected$Super_Bad = collected$Spearman > quantile(collected$Spearman, probs = 2/3) & 
                      collected$Gleason_Morisita > median(collected$Gleason_Morisita)

# Calculate fit for tertile split in Spearman instead
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,")~ Super_Bad")), data = collected)

# Plot most heterogeneous section
p = ggsurvplot(fit,
               pval = TRUE, conf.int = FALSE,
               legend = "none",
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = c(1,1), # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               tables.y.text = FALSE,
               ylab = paste0("Fraction ",out_var," Free"),
               xlab = "Time (Years)",
               xscale = 12,
               break.x.by = 12*4,
               tables.theme = theme(axis.title.x=element_blank(),
                                    axis.text.x=element_blank(),
                                    axis.ticks.x=element_blank(),
                                    plot.title = element_blank()),
               palette = c("#929591", "#e50000"))
ggsave(plot = print(p), "results/9_outcome_analysis/Spearman_tert_GM_median_survival.png", width = 4, height = 4)

# Read in the geometric mean
collected$Geometric_Mean = scaleQuant(sqrt(collected$Spearman*collected$Gleason_Morisita))

# Make categorical
collected$High_spearman = collected$Spearman > quantile(collected$Spearman, probs = 2/3)
collected$High_GM       = collected$Gleason_Morisita > median(collected$Gleason_Morisita)

# What do we want to look at?
covariates = c("Geometric_Mean")

covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", covariates)

# Other variables
collected$High_PSA               = collected$Presenting_PSA > 20
collected$Is_T3                  = collected$T_stage == "T3"
collected$N_stage                = collected$N_stage == "N1-3"
collected$Gleason_grade_group    = as.factor(collected$Gleason_Grade_Rescore)

# To make 5 the reference
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 5:2)

# Now, using the same covariates, perform a multivariate test
res.cox = coxph(as.formula(paste0("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected)

# For plotting
collected$Gleason_grade_group     = factor(collected$Gleason_grade_group, levels = 2:5)

# Collect the results for a table 
sum.res.cox = summary(res.cox)

# Format
res.table = signif(cbind(sum.res.cox[["conf.int"]][,c(1,3,4)], sum.res.cox[["coefficients"]][,"Pr(>|z|)"]), digits = 3)

# To make it readable
colnames(res.table) = c("HR", "Lower", "Upper", "p")
rownames(res.table) = c("Samples", "PSA_20+", "Gleason_4", "Gleason_3", "Gleason_2", "T3+", "N1+", "Joint Diversity")

# Change R format
res.table = as.data.frame(res.table)

# To be read!
res.table$Reference = c("0", "PSA_less_than_20", 
                        "Gleason_5", "Gleason_5", "Gleason_5", 
                        "Not_T3", "N0", "0")

# Write it out
write.csv(res.table, file = "results/9_outcome_analysis/SpearGMGM_CPH_model.csv", quote = F)

# Use survminer to plot the results
p = ggforest(res.cox)
ggsave(plot = print(p), "results/9_outcome_analysis/SpearGMGM_CPH_model.png", width = 8, height = 7)
ggsave(plot = p, "results/9_outcome_analysis/SpearGMGM_CPH_model.pdf", width = 8, height = 7)
