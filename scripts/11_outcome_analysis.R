# This is a script to perform outcome analysis

library(survival)
library(survminer)
library(pROC)
library(ggplot2)
library(cowplot)

# Source general scripts
source("scripts/00_general_functions.R")

##### The parameters #####

# Outcome
out_var = "Recurrence"
out_tim = paste0(out_var,"_Months")

# Spearman or geometric mean of spearman and gleason morisita (only for testing)
spearman_choice = c("Spearman", "SpearGMGM")[1]

# Parameters for clinical variables
cap_time      = Inf
psa_threshold = 20

# p value threshold for being included in the multivariate
cont_uni_pv   = 0.1
cate_uni_pv   = 0.1

# Which threshold to use?
thresh_mode = c("roc_threshold", "median", "tertile", "upper_tert")[2]

output_report = F
plot_width    = 12

group_palette = c("#0343df", "#f97306", "#e50000")

#################################################################################################################

# Read in the outcome data
outcome = read.csv("refs/FORECAST_outcome_data_v2.csv", stringsAsFactors = F)

# Read in the data characteristics 
characteristics = read.csv("refs/FORECAST_patient_characteristics.csv", stringsAsFactors = F)

# Read in Gleason data
digipath_gleason = read.csv("refs/GleasonGradeGroupPatient.csv", 
                            header = F, stringsAsFactors = F)
colnames(digipath_gleason) = c("Patient", "Digital_Gleason")
digipath_gleason$Digital_Gleason = floor(digipath_gleason$Digital_Gleason)

# Read in rescoring pathologist
gleason_rescore = read.csv("refs/GleasonRescore_patients.csv", 
                           header = T, stringsAsFactors = F)

# Read in the metrics collected
metrics = read.table("results/8_collect_metrics/FORECAST_genomic_metrics.txt", sep = "\t", header = T)

# Read in GM data too
gleason_morisita = read.csv("refs/GleasonMorisita_Patients.csv", 
                            header = F, stringsAsFactors = F)

# Name GM columns
colnames(gleason_morisita) = c("Patient", "Gleason_Morisita")

# Make columns more readable
colnames(outcome) = c("Patient", "Recurrence", "Recurrence_Months", "Death", "Death_Months", "Mets", 
                      "Mets_Months", "Wide_Mets", "Wide_Met_Months")

# Rename the characteristics column
colnames(characteristics)[1] = "Patient"

# Cap?
outcome[,out_var]        = ifelse(outcome[,out_tim] > cap_time, yes = 0, no = outcome[,out_var])
outcome[,out_tim]        = ifelse(outcome[,out_tim] > cap_time, yes = cap_time, no = outcome[,out_tim])

# Add IM to beginning 
outcome$Patient         = ifelse(grepl("DNT", outcome$Patient), 
                                 yes = outcome$Patient, no = paste0("IM",outcome$Patient))
characteristics$Patient = ifelse(grepl("DNT", characteristics$Patient), 
                                 yes = characteristics$Patient, no = paste0("IM",characteristics$Patient))

# Merge by patient id
collected = merge(outcome, characteristics, by = "Patient")
collected = merge(collected, metrics, by = "Patient")
collected = merge(collected, gleason_morisita, by = "Patient")
collected = merge(collected, digipath_gleason, by = "Patient")
collected = merge(collected, gleason_rescore[,c("Patient", "Gleason_Grade_Rescore")], by = "Patient", all.x = T)

# Remove delineate samples
collected = collected[grepl("IM", collected$Patient),]

# Remove anything that doesn't have enough samples
collected = collected[collected$Samples > 2,]

# Make into categories
collected$High_PSA               = collected$Presenting_PSA > psa_threshold
collected$Is_T3                  = collected$T_stage == "T3"
collected$N_stage                = collected$N_stage == "N1-3"
collected$Gleason_grade_group    = as.factor(collected$Gleason_Grade_Rescore)

# Add Geometric Mean of Spearman and GM
collected$SpearGMGM = sqrt(collected$Spearman*collected$Gleason_Morisita)

# Gleason grade group and PGA
ggplot(collected, aes(x = Gleason_grade_group, y = mPGA)) + xlab("Gleason grade group") + geom_violin() + 
  geom_boxplot(width=0.1) + theme_cowplot()

if(output_report) {
  # Save in a mega pdf because there are lots of comparison
  pdf(paste0("results/11_outcome_analysis/",out_var,
             "_timelimit=",cap_time,"_thresholdtype=",thresh_mode,"_outcome_analysis_report.pdf"), 
      width = 7, height = 7)
}

if(thresh_mode=="median") {
  
  collected$High_mPGA             = collected$mPGA > median(collected$mPGA)
  collected$High_sdPGA            = collected$sdPGA > median(collected$sdPGA)
  collected$High_max_PGA          = collected$max_PGA > median(collected$max_PGA)
  collected$High_subclonal_events = collected$Number_subclonal_events > median(collected$Number_subclonal_events)
  collected$High_SpearGMGM        = collected$SpearGMGM > median(collected$SpearGMGM)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_mPGA             = tertile_split(collected$mPGA)
  collected$High_sdPGA            = tertile_split(collected$sdPGA)
  collected$High_max_PGA          = tertile_split(collected$max_PGA)
  collected$High_subclonal_events = tertile_split(collected$Number_subclonal_events)
  collected$High_SpearGMGM        = tertile_split(collected$SpearGMGM)
  
}

if(thresh_mode=="upper_tert") {
  
  collected$High_mPGA             = collected$mPGA > quantile(collected$mPGA, probs = 2/3)
  collected$High_sdPGA            = collected$sdPGA > quantile(collected$sdPGA, probs = 2/3)
  collected$High_max_PGA          = collected$max_PGA > quantile(collected$max_PGA, probs = 2/3)
  collected$High_subclonal_events = collected$Number_subclonal_events > quantile(collected$Number_subclonal_events, probs = 2/3)
  collected$High_SpearGMGM        = collected$SpearGMGM > quantile(collected$SpearGMGM, probs = 2/3)
  
}

# Calculate fit for mean PGA
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_mPGA")), data = collected)

# Plot mean PGA result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate fit for sd PGA
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_sdPGA")), data = collected)

# Plot sd PGA result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate fit for max PGA
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_max_PGA")), data = collected)

# Plot max PGA result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate fit for number of subclonal events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_subclonal_events")), data = collected)

# Plot the fit for number of subclonal events
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate fit for TP53, removing the NAs
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ TP53")), 
              data = colNArm(collected, sub = "TP53"))

# Plot the TP53 result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate fit for ATM, removing the NAs
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ ATM")), 
              data = colNArm(collected, sub = "TP53"))

# Plot the ATM result
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate whether a "DNA damage" mutation has happened
collected$DDM_status = collected$DNA_damage_mutation!="Wild-type"

# Calculate the fit for these DDMs
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ DDM_status")), 
              data = colNArm(collected, sub = "DDM_status"))

# Plot the fit for these mutations
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a subclonal mutation or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ Subclonal_Mut_Status")), 
              data = colNArm(collected, sub = "Subclonal_Mut_Status"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a MYCN amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYCN")), 
              data = colNArm(collected, sub = "MYCN"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a MDM2 amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MDM2")), 
              data = colNArm(collected, sub = "MDM2"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a MYC amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYC")), 
              data = colNArm(collected, sub = "MYC"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a FGFR1 amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ FGFR1")), 
              data = colNArm(collected, sub = "FGFR1"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a key amplification or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ key_amp")), 
              data = colNArm(collected, sub = "key_amp"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a MYCN or MDM2 amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYCN_MDM2")), 
              data = colNArm(collected, sub = "key_amp"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

# Calculate the fit for the patient's having a MYC or FGFR1 amp or not
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ MYC_FGFR1")), 
              data = colNArm(collected, sub = "key_amp"))

# Plot that fit too
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(thresh_mode=="median") {
  
  collected$High_Gleason = collected$Cont_Gleason > median(collected$Cont_Gleason, na.rm = T)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_Gleason = tertile_split(collected$Cont_Gleason)
  
}

if(thresh_mode=="upper_tert") {
  
  collected$High_Gleason = collected$Cont_Gleason > quantile(collected$Cont_Gleason, probs = 2/3, na.rm = T)
  
}

# Make a fit object
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_Gleason")), 
              data = colNArm(collected, sub = "High_Gleason"))

# Plot results for gleason score
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(thresh_mode=="median") {
  
  collected$High_Tot_Events = collected$Total_Events > median(collected$Total_Events)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_Tot_Events = tertile_split(collected$Total_Events)
  
}

if(thresh_mode=="upper_tert") {
  
  collected$High_Tot_Events = collected$Total_Events > quantile(collected$Total_Events, probs = 2/3)
  
}

# Make a fit object for total events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_Tot_Events")), data = collected)

# Plot the results of the total event threshold
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(thresh_mode=="median") {
  
  collected$High_subclonality = collected$Subclonality > median(collected$Subclonality)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_subclonality = tertile_split(collected$Subclonality)
  
}

if(thresh_mode=="upper_tert") {
  
  collected$High_subclonality = collected$Subclonality > quantile(collected$Subclonality, probs = 2/3)
  
}

# Make a fit object for the % of subclonal events
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_subclonality")), data = collected)

# Plot that result out 
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(thresh_mode=="median") {
  
  collected$High_lossiness = collected$Lossness > median(collected$Lossness)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_lossiness = tertile_split(collected$Lossness)
  
}

if(thresh_mode=="upper_tert") {
  
  collected$High_lossiness = collected$Lossness > quantile(collected$Lossness, probs = 2/3)
  
}

# Make a fit object for lossness
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_lossiness")), data = collected)

# Plot that result out 
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(thresh_mode=="median") {
  
  collected$High_spearman = collected$Spearman > median(collected$Spearman)
  
}

if(thresh_mode=="tertile") {
  
  collected$High_spearman = tertile_split(collected$Spearman)
  
}

if(thresh_mode=="upper_tert") {
  
  up_tert_threshold = quantile(collected$Spearman, probs = 2/3)
  
  collected$High_spearman = collected$Spearman > up_tert_threshold
  ggplot(collected, aes(x = Spearman)) + geom_histogram(color="black", fill="white", binwidth = 0.01) + 
    geom_vline(xintercept = up_tert_threshold, lty = "dotted") + 
    scale_y_continuous(name = "Count", limits = c(0,10), breaks = seq(0,10,by = 2)) + 
    ggtitle(paste0("Spearman values for IMRT patients used for outcome analysis\n(upper tertile split = ",
                   signif(up_tert_threshold, digits = 3),")")) + theme_cowplot() + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave("results/11_outcome_analysis/Spearman_distribution_tertile_split.png", width = 7.5, height = 7.5)
  
}

# Make a fit object for the Spearman metric
fit = survfit(as.formula(paste0("Surv(",out_tim,", ",out_var,") ~ High_spearman")), data = collected)

# Plot that result out 
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = group_palette)

if(output_report) {dev.off()}

#################################################################################################################
######################################## Cox Proportional Hazards Model #########################################
#################################################################################################################

# Half ploidy to make it tetraploid or not as a measure
collected$Ploidy = collected$Ploidy / 2

# Transform some variables
collected$mPGA                    = log(collected$mPGA)
collected$sdPGA                   = log(collected$sdPGA)
collected$max_PGA                 = log(collected$max_PGA)
collected$Subclonality            = exp(collected$Subclonality)
collected$Lossness                = log(collected$Lossness)
collected$Number_subclonal_events = log(collected$Number_subclonal_events)
collected$Total_Events            = log(collected$Total_Events)

# Bug fix
if(spearman_choice=="Spearman") {
  threshold_spearman_choice = "High_spearman"
} else {threshold_spearman_choice = paste0("High_",spearman_choice)}

collected$TP53_category       = factor(collected$TP53_category, 
                                       levels = c("Wild-type", "Subclonal", "Clonal"))
collected$DNA_damage_mutation = factor(collected$DNA_damage_mutation, 
                                       levels = c("Wild-type", "Subclonal", "Clonal"))
collected$ATM_category        = factor(collected$ATM_category, 
                                       levels = c("Wild-type", "Subclonal", "Clonal"))

# What do we want to investigate?
covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", "mPGA", "sdPGA", "max_PGA", "L2RSS", 
               "Lossness", "Total_Events", "Subclonality", "Number_subclonal_events", 
               "Subclonal_Mut_Status", "Ploidy", "DNA_damage_mutation", "key_amp", "MYCN_MDM2", "MYC_FGFR1", spearman_choice)

# Get the formulas for each univariate test
univ_formulas = sapply(covariates,
                       function(x) as.formula(paste0("Surv(",out_tim,", ",out_var,")~", x)))

# Perform univariate CPH model
univ_models = lapply(univ_formulas, function(x){coxph(x, data = collected)})

# Extract data 
univ_results = lapply(univ_models, function(x) {
  
  x = summary(x)
  p.value = signif(x$wald["pvalue"], digits=2)
  wald.test = signif(x$wald["test"], digits=2)
  beta = signif(x$coef[1], digits=2);#coeficient beta
  HR = signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  HR = paste0(HR, " (", 
              HR.confint.lower, "-", HR.confint.upper, ")")
  res = c(beta, HR, wald.test, p.value)
  names(res) = c("beta", "HR (95% CI for HR)", "wald.test", 
                 "p.value")
  return(res)

})

# Which qualify for the next round?
covariates = covariates[unlist(lapply(univ_models, function(i) any(summary(i)$coefficients[,5] < cont_uni_pv)))]

covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", covariates)

# Now, using the same covariates, perform a multivariate test
res.cox = coxph(as.formula(paste("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected)

# Use survminer to plot the results
ggforest(res.cox)
ggsave(paste0("results/11_outcome_analysis/",out_var,
              "_timelimit=",cap_time,
              "_thresholdtype=",thresh_mode,"_outcome_analysis_Cox_proportional_hazards_model_continuous.pdf"), 
       width = plot_width, height = 7)

# What do we want to investigate?
covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", 
               "High_mPGA", 
               "High_sdPGA",
               "High_max_PGA", 
               "High_subclonal_events", 
               "High_Tot_Events", 
               "High_subclonality",
               "High_lossiness", 
               threshold_spearman_choice,
               "Subclonal_Mut_Status",
               "Ploidy", 
               "DNA_damage_mutation",
               "key_amp", "MYCN_MDM2", "MYC_FGFR1")

# Get the formulas for each univariate test
univ_formulas = sapply(covariates,
                       function(x) as.formula(paste0("Surv(",out_tim,", ",out_var,")~", x)))

# Perform univariate CPH model
univ_models = lapply(univ_formulas, function(x){coxph(x, data = collected)})

# Extract data 
univ_results = lapply(univ_models, function(x) {
  
  x = summary(x)
  p.value = signif(x$wald["pvalue"], digits=2)
  wald.test = signif(x$wald["test"], digits=2)
  beta = signif(x$coef[1], digits=2);#coeficient beta
  HR = signif(x$coef[2], digits=2);#exp(beta)
  HR.confint.lower = signif(x$conf.int[,"lower .95"], 2)
  HR.confint.upper = signif(x$conf.int[,"upper .95"],2)
  HR = paste0(HR, " (", 
              HR.confint.lower, "-", HR.confint.upper, ")")
  res = c(beta, HR, wald.test, p.value)
  names(res) = c("beta", "HR (95% CI for HR)", "wald.test", 
                 "p.value")
  return(res)
  
})

# Which ones qualify?
covariates = covariates[unlist(lapply(univ_models, function(i) any(summary(i)$coefficients[,5] < cate_uni_pv)))]

covariates = c("Samples", "High_PSA", "Gleason_grade_group", "Is_T3", "N_stage", covariates)

# Now, using the same covariates, perform a multivariate test
res.cox = coxph(as.formula(paste0("Surv(",out_tim,", ",out_var,")~", paste0(covariates, collapse = " + "))), 
                data = collected)

# Use survminer to plot the results
ggforest(res.cox)
ggsave(paste0("results/11_outcome_analysis/",
              out_var,"_timelimit=",cap_time,
              "_thresholdtype=",thresh_mode,"_outcome_analysis_Cox_proportional_hazards_model_categorical.pdf"), 
       width = plot_width, height = 7)
