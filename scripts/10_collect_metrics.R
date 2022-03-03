# Collect our metrics
library(plyr)
library(ggplot2)
library(cowplot)
library(GGally)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
theme_set(theme_cowplot())

# Functions
source("scripts/00_general_functions.R")

# Minimum PGA
min_pga = 0.01

# Minimum number of samples
min_samples = 3

# Targeted data repo
targeted_repo = "target_repo/"

# Read info on tes patients
tes_patients = read.table(paste0(targeted_repo,"results/mutation_calling/Patients_with_normal_attempted_mut.txt"))
insuff_sampl = read.table(paste0(targeted_repo,"results/mutation_calling/Patients_excluded_TES_too_few_samples.txt"))
exclude_sm   = read.table(paste0(targeted_repo,"results/mutation_calling/Patients_excluded_TES_exclude.txt"))
insuff_sampl = rbind(insuff_sampl, exclude_sm)

# Read in the pipeline metrics
metrics    = read.table("results/2_metrics/metric_assessment_table.txt", header = T, stringsAsFactors = F)

# Pre-filtering count
pre_metrics = metrics[!(metrics$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", metrics$Sample) | metrics$Sample=="IM1017_B2_DW1"),]
pre_metrics = ddply(pre_metrics, .(Patient), 
                    summarize, 
                    Prefilter_samples = length(PGA))

# Remove bad samples
metrics = metrics[!(metrics$PGA < min_pga | metrics$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", metrics$Sample) | metrics$Sample=="IM1017_B2_DW1"),]

# Summarise
pp_metrics = ddply(metrics, .(Patient), 
                   summarize, 
                   Samples = length(PGA), 
                   mPGA = mean(PGA), 
                   sdPGA = sd(PGA),
                   max_PGA = max(PGA))

# Add that prefilter data to the file
pp_metrics = merge(pre_metrics, pp_metrics)

# Read in the l2rss data
l2rss = readRDS("results/5_boxplot_log2ratio_dist/L2RSS_scores_per_patient.rds")
l2rss = list_summary_df(l2rss, value_name = "L2RSS", summary = "median")

# Read in the spearman data
spear = readRDS("results/5_boxplot_log2ratio_dist/Spearman_rho_per_patient.rds")
spear = list_summary_df(spear, value_name = "Spearman", summary = "mean")

# Read in lossness data
lossness = readRDS("results/6_lossness/FORECAST_collapsed_segment_lossness.rds")
lossness = list_summary_df(lossness, value_name = "Lossness", summary = "mean")

# Read in lossness data
medicc = read.table("results/7_medicc_analysis/medicc_analysis_summary.txt", header = T, sep = "\t")
medicc = medicc[,c("Patient", "Total_Events", "Subclonality")]
medicc$Number_subclonal_events = medicc$Total_Events * medicc$Subclonality

# Read in Dropbox
morisita = read.csv("refs/Morisita_Cancer_Regions.csv", 
                    stringsAsFactors = F)

# Edit patient/sample names
morisita$Sample = convertPathSampleName(morisita$SampleID)

# Read in the lymphocyte stuff
lymphocyte = read.csv("refs/ITLR_Cancer_Regions.csv", 
                      stringsAsFactors = F)

# Get the correct name
lymphocyte$Sample = convertPathSampleName(lymphocyte$X)

# Read in cell class stuff
cell_class = read.csv("refs/CellClassificationSummary_Regions.csv",
                      stringsAsFactors = F)

# Change name
cell_class$Sample = convertPathSampleName(cell_class$X)

# Merge the data
cc_collected = merge(metrics[,c("Sample", "Patient")], morisita, by = "Sample")
cc_collected = merge(cc_collected, lymphocyte, by = "Sample")
cc_collected = merge(cc_collected, cell_class, by = "Sample")

# Summarise
cc_metrics = ddply(cc_collected, .(Patient), 
                   summarize,
                   mMorisita = mean(Morisita),
                   mATLR = mean(ATLR),
                   mITLR = mean(ITLR),
                   mICs = mean(Immune.Clustering))

# Read in gleason
gleason = read.csv("refs/GleasonSummaryClean.csv", 
                   stringsAsFactors = F, header = T)
gleason$Patient = unlist(lapply(strsplit(gleason$X, split = ""), function(i) paste0(i[1:6], collapse = "")))

# Read in the new per region data
gleason_regions = read.csv("refs/GleasonSummary_Regions.csv", 
                           stringsAsFactors = F, header = T)
gleason_regions$Patient   = unlist(lapply(strsplit(gleason_regions$X, split = ""), function(i) paste0(i[1:6], collapse = "")))
gleason_regions$Region   = unlist(lapply(strsplit(gleason_regions$X, split = ""), function(i) paste0(i[-c(1:6, length(i))], collapse = "")))
gleason_regions$Subregion = unlist(lapply(strsplit(gleason_regions$X, split = ""), function(i) paste0(i[length(i)], collapse = "")))

# Get conversion reference table
look_up_table_slides = read.csv("refs/Forecast ID Reference.csv", 
                                stringsAsFactors = F, header = T)

# Filter for those we have sequencing data
sequenced_slides = look_up_table_slides[look_up_table_slides$Region!="",]
sequenced_slides$Region = gsub(" ", "", sequenced_slides$Region)
sequenced_slides$X = paste0(sequenced_slides$Patient,sequenced_slides$Slide)

# Merge gleason and slide data
gleason_slide = merge(gleason, sequenced_slides[,2:5], by = "X")
gleason_slide$Block = NULL
gleason_slide$Slide = NULL
gleason_slide$Subregion = ""

# Now for both do a full name
gleason_regions$Sample = paste0(gleason_regions$Patient,"_",gleason_regions$Region,gleason_regions$Subregion,"_DW1")
gleason_slide$Sample = paste0(gleason_slide$Patient,"_",gleason_slide$Region,gleason_slide$Subregion,"_DW1")

# Get only sequenced regions
gleason_all = rbind(gleason_slide, gleason_regions)

# Now replace what we called gleason with the new df
gleason = gleason_all[gleason_all$Sample %in% metrics$Sample,]

# Collapse
gleason_patient = ddply(gleason, .(Patient), summarize,
                        sumGS3 = sum(Gleason.3),
                        sumGS4 = sum(Gleason.4),
                        sumGS5 = sum(Gleason.5))

gleason_patient$GS_count = apply(gleason_patient, 1, function(i) which.max(i[3:5])+2)

# Add a continuous gleason score
gleason_patient$Cont_Gleason = apply(gleason_patient, 1, function(i) {
  
  sum(as.numeric(i[c("sumGS3", "sumGS4", "sumGS5")])*c(3:5)) / sum(as.numeric(i[c("sumGS3", "sumGS4", "sumGS5")]))
  
})

# Get ploidy 
ploidy = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                    header = T, sep = "\t", stringsAsFactors = F)
colnames(ploidy)[1] = "Patient"

# Round it off
ploidy$Ploidy = round(ploidy$PsiT)

# Get subclonal mutation status
subclonal_mut_patient = read.table(paste0(targeted_repo,"results/mutation_calling/Subclonal_driver_status_summary.txt"), 
                                   sep = "\t", header = T, stringsAsFactors = F)

# Remove those where we don't have enough data to say if it's subclonal
subclonal_mut_patient$Subclonal_Mut_Status[subclonal_mut_patient$Patient %in% insuff_sampl$V1] = NA

# Read it in
amps_detection = readRDS("results/3_amplification_detection/subclonality_assessment.rds")

# Detect different amps
amps_collected = data.frame(Patient = names(amps_detection),
                            MYCN = unlist(lapply(amps_detection, function(i) "MYCN" %in% names(i))),
                            MDM2 = unlist(lapply(amps_detection, function(i) "MDM2" %in% names(i))),
                            MYC  = unlist(lapply(amps_detection, function(i) "MYC" %in% names(i))),
                            FGFR1  = unlist(lapply(amps_detection, function(i) "FGFR1" %in% names(i))))
rownames(amps_collected) = NULL

# Add flag for any amplification
amps_collected$key_amp   = apply(amps_collected[,2:ncol(amps_collected)], 1, any)
amps_collected$MYCN_MDM2 = apply(amps_collected[,c("MYCN", "MDM2")], 1, any)
amps_collected$MYC_FGFR1 = apply(amps_collected[,c("MYC", "FGFR1")], 1, any)

# Do the merge
collected = merge(pp_metrics, l2rss, by = "Patient")
collected = merge(collected, lossness, by = "Patient")
collected = merge(collected, spear, by = "Patient")
collected = merge(collected, gleason_patient, by = "Patient", all.x = T)
collected = merge(collected, medicc, by = "Patient")
collected = merge(collected, cc_metrics, by = "Patient", all.x = T)
collected = merge(collected, subclonal_mut_patient, by = "Patient")
collected = merge(collected, amps_collected, by = "Patient")

# Transform spearman
collected$Spearman = 1 - collected$Spearman

# Make a plot
p = ggpairs(collected[,!colnames(collected) %in% c("medianGS_ps", "Patient", "sumGS3", 
                                                   "sumGS4", "sumGS5", "Subclonal_Mut_Status",
                                                   "MYCN", "MDM2", "MYC", "FGFR1", 
                                                   "key_amp", "MYCN_MDM2", "MYC_FGFR1")])
ggsave(plot = p, filename = "results/8_collect_metrics/FORECAST_cohort_metrics.pdf", width = 13, height = 10)

# Normalise and make matrix
mat = collected[,3:ncol(collected)]

# Read in the mutations
muts = readRDS(paste0(targeted_repo,"results/mutation_calling/subclonality_assessment.rds"))

# Annotate
right_annotation = rowAnnotation(Gleason = collected$GS_count,
                                 TP53 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$TP53), 
                                                       attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1),
                                 SPOP = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$SPOP), 
                                                       attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1),
                                 CTNNB1 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$CTNNB1), 
                                                         attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1),
                                 AKT1 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$AKT1), 
                                                       attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1),
                                 col = list(Gleason = c("3" = "#15b01a",
                                                        "4" = "#0343df",
                                                        "5" = "#f97306"),
                                            TP53 = c("FALSE" = "#6b8ba4", "TRUE" = "#dbb40c"),
                                            SPOP = c("FALSE" = "#6b8ba4", "TRUE" = "#dbb40c"),
                                            CTNNB1 = c("FALSE" = "#6b8ba4", "TRUE" = "#dbb40c"),
                                            AKT1 = c("FALSE" = "#6b8ba4", "TRUE" = "#dbb40c")))

mat = apply(mat, 2, function(i) {(i - min(i, na.rm = T)) / (max(i, na.rm = T) - min(i, na.rm = T))})
mat = mat[,!colnames(mat) %in% c("GS_count", "medianGS_ps", "sumGS3", "sumGS4", "sumGS5", "Cont_Gleason")]
Heatmap(mat, col = colorRamp2(c(0, 1), c("white", "red")), right_annotation = right_annotation, 
        name = "Normalised\nValue")

collected = merge(collected, ploidy, by = "Patient")

#################################################################################################################
################ Add in all the things to the collected object we want ##########################################
#################################################################################################################

# What's the TP53 status?
collected$TP53 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$TP53), 
                                attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# What's the BRCA1/2 status?
collected$BRCA = applyMutStatus(query_list = collected$Patient, mutated_cases = unique(c(names(muts$BRCA1), 
                                                                                         names(muts$BRCA2))), 
                                attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# What's the BRCA1/2 clonal status?
collected$BRCA_clonal = applyMutStatus(query_list = collected$Patient, 
                                       mutated_cases = unique(c(names(which(unlist(lapply(muts$BRCA1, function(i) any(!i))))),
                                                                names(which(unlist(lapply(muts$BRCA2, function(i) any(!i))))))),
                                       attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# BRCA categories
collected$BRCA_category = ifelse(collected$BRCA_clonal, "Clonal", ifelse(collected$BRCA & !collected$BRCA_clonal,
                                                                         yes = "Subclonal", no = "Wild-type"))

# What's the ATM clonal status?
collected$ATM_clonal = applyMutStatus(query_list = collected$Patient, 
                                      mutated_cases = names(which(unlist(lapply(muts$ATM, function(i) any(!i))))), 
                                      attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# What's the ATM status?
collected$ATM = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$ATM), 
                               attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# ATM categories
collected$ATM_category = ifelse(collected$ATM_clonal, "Clonal", ifelse(collected$ATM & !collected$ATM_clonal, 
                                                                       yes = "Subclonal", no = "Wild-type"))

# What's the PALB2 status?
collected$PALB2 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$PALB2), 
                                 attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# What's the PALB2 clonal status?
collected$PALB2_clonal = applyMutStatus(query_list = collected$Patient, 
                                        mutated_cases = names(which(unlist(lapply(muts$PALB2, function(i) any(!i))))),
                                        attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# PALB2 categories
collected$PALB2_category = ifelse(collected$PALB2_clonal, "Clonal", ifelse(collected$PALB2 & !collected$PALB2_clonal,
                                                                           yes = "Subclonal", no = "Wild-type"))

# What's the CDK12 status?
collected$CDK12 = applyMutStatus(query_list = collected$Patient, mutated_cases = names(muts$CDK12), 
                                 attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# What's the CDK12 clonal status?
collected$CDK12_clonal = applyMutStatus(query_list = collected$Patient, 
                                        mutated_cases = names(which(unlist(lapply(muts$CDK12, function(i) any(!i))))),
                                        attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# CDK12 categories
collected$CDK12_category = ifelse(collected$CDK12_clonal, "Clonal", ifelse(collected$CDK12 & !collected$CDK12_clonal,
                                                                           yes = "Subclonal", no = "Wild-type"))

# What's the TP53 status?
collected$TP53_clonal = applyMutStatus(query_list = collected$Patient, 
                                       mutated_cases = names(which(unlist(lapply(muts$TP53, function(i) any(!i))))), 
                                       attempted_seq = tes_patients$V1, bad_cases = insuff_sampl$V1)

# TP53 categories
collected$TP53_category = ifelse(collected$TP53_clonal, "Clonal", ifelse(collected$TP53 & !collected$TP53_clonal,
                                                                         yes = "Subclonal", no = "Wild-type"))

# DNA damage mutation status - if any is clonal, it's classified as clonal
collected$DNA_damage_mutation = apply(collected[,c("TP53_category", "ATM_category", 
                                                   "BRCA_category", "CDK12_category", 
                                                   "PALB2_category")], 1, function(i) {

  if(!all(is.na(i))) {
  
    res = "Subclonal"
    
    if(any(i == "Clonal")) {res = "Clonal"}
    if(all(i == "Wild-type")) {res = "Wild-type"}
  
  } else {res = NA}
  
  return(res)

})

# Write out all that lovely data
write.table(collected, file = "results/8_collect_metrics/FORECAST_genomic_metrics.txt", 
            quote = F, sep = "\t", row.names = F)

#################################################################################################################
################ Now do all the plots with the sample filter ####################################################
#################################################################################################################

# Subset now by the minimum number of samples required
collected = collected[collected$Samples >= min_samples,]

# P Values
p_value = signif(t.test(mPGA ~ TP53, colNArm(collected, sub = "TP53"))$p.value, digits = 3)
p_value_line_height = max(collected$mPGA)*1.1

ggplot(colNArm(collected, sub = "TP53"), aes(x = TP53, y = mPGA)) + geom_boxplot() + 
  ggtitle(paste0("TP53 mutation status and mPGA")) +
  xlab("TP53 mutated?") + 
  ylab("Mean PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_PGA_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(max_PGA ~ TP53, colNArm(collected, sub = "TP53"))$p.value, digits = 3)
p_value_line_height = max(collected$max_PGA)*1.1

ggplot(colNArm(collected, sub = "TP53"), aes(x = TP53, y = max_PGA)) + geom_boxplot() + 
  ggtitle(paste0("TP53 mutation status and max PGA")) +
  xlab("TP53 mutated?") + 
  ylab("Max PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_max_PGA_correlation.pdf", height = 6, width = 7)

ggplot(colNArm(collected, sub = "TP53"), aes(x = TP53, y = PsiT)) + geom_boxplot() + 
  ggtitle(paste0("TP53 mutation status and PsiT")) +
  xlab("TP53 mutated?") + 
  ylab("Fitted Ploidy") + 
  ylim(1.5,4.5) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_PsiT_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(Lossness ~ BRCA, colNArm(collected, sub = "BRCA"))$p.value, digits = 3)
p_value_line_height = max(collected$Lossness)*1.1

ggplot(colNArm(collected, sub = "BRCA"), aes(x = BRCA, y = Lossness)) + geom_boxplot() + 
  ggtitle(paste0("BRCA1/2 mutation status and Lossiness")) +
  xlab("BRCA1/2 mutated?") + 
  ylab("Lossiness") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/BRCA_Lossness_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(mPGA ~ ATM_clonal, colNArm(collected, sub = "ATM_clonal"))$p.value, digits = 3)
p_value_line_height = max(collected$mPGA)*1.1

ggplot(colNArm(collected, sub = "ATM_clonal"), aes(x = ATM_clonal, y = mPGA)) + geom_boxplot() + 
  ggtitle(paste0("Clonal ATM mutation status and mPGA")) +
  xlab("Clonal ATM mutation?") + 
  ylab("mPGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/ATM_clonal_mPGA_correlation.pdf", height = 6, width = 7)

# Compare PGA and Events
ggplot(collected, aes(x = mPGA, y = Total_Events)) + geom_point() + 
  ggtitle(paste0("Mean PGA v Total Events")) +
  xlab("Mean PGA") + 
  ylab("Total Events") + 
  ylim(0,300) + 
  xlim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/PGA_Events_correlation.pdf", height = 6, width = 7)

# p values
p_value = signif(summary(lm(Subclonality ~ mPGA, data = collected))$coefficients[2,4], digits = 3)

# Compare PGA and Events
ggplot(collected, aes(x = mPGA, y = Subclonality)) + geom_point() + 
  ggtitle(paste0("Mean PGA v Subclonality (p=",p_value,")")) +
  xlab("Mean PGA") + 
  ylab("Subclonality") + 
  ylim(0,1) + 
  xlim(0,1) + theme_cowplot() + geom_smooth(method = "lm") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/PGA_Subclonality_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(Total_Events ~ TP53, colNArm(collected, sub = "TP53"))$p.value, digits = 4)
p_value_line_height = max(colNArm(collected, sub = "TP53")$Total_Events)*1.1

# Let's have a look at events
ggplot(colNArm(collected, sub = "TP53"), aes(x = TP53, y = Total_Events)) + geom_boxplot() + 
  ggtitle(paste0("TP53 mutation status and tree events")) +
  xlab("TP53 mutation?") + 
  ylab("Total Events") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height * 0.99, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height * 1.05, size = 5) +
  ylim(0,350) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_Events_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(summary(lm(mPGA ~ Cont_Gleason, data = collected))$coefficients[2,4], digits = 3)

ggplot(collected, aes(x = Cont_Gleason, y = mPGA)) + geom_point() + 
  ggtitle(paste0("Mean PGA v Gleason (p=",p_value,")")) +
  ylab("Mean PGA") + 
  xlab("Continuous Gleason") + 
  xlim(3,5) + 
  ylim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_mPGA_correlation.pdf", height = 6, width = 7)

ggplot(collected, aes(x = Cont_Gleason, y = mPGA)) + geom_point() +
  ylab("Mean PGA") +
  xlab("Continuous Gleason") +
  xlim(3,5) +
  ylim(0,0.6) +
  geom_smooth(method = "lm", col = "#0165fc", se = T) +
  annotate(geom="text", x=3.5, y=0.45, label=paste0("p=",format(p_value, scientific = FALSE)), size = 5) +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_mPGA_correlation.png", height = 3, width = 4)


ggplot(collected, aes(x = mPGA, y = mMorisita)) + geom_point() + 
  #ggtitle(paste0("Mean PGA v Gleason (p=",p_value,")")) +
  xlab("Mean PGA") + 
  ylab("Mean Morisita Index") + 
  xlim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/mMorisita_mPGA_correlation.pdf", height = 6, width = 7)


ggplot(collected, aes(x = mPGA, y = mATLR)) + geom_point() + 
  #ggtitle(paste0("Mean PGA v Gleason (p=",p_value,")")) +
  xlab("Mean PGA") + 
  ylab("Mean Adjacent Tumour Lymphocyte Ratio") + 
  xlim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/mATLR_mPGA_correlation.pdf", height = 6, width = 7)


ggplot(collected, aes(x = mPGA, y = mITLR)) + geom_point() + 
  #ggtitle(paste0("Mean PGA v Gleason (p=",p_value,")")) +
  xlab("Mean PGA") + 
  ylab("Mean Inta-Tumour Lymphocyte Ratio") + 
  xlim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/mITLR_mPGA_correlation.pdf", height = 6, width = 7)


ggplot(collected, aes(x = mPGA, y = mICs)) + geom_point() + 
  #ggtitle(paste0("Mean PGA v Gleason (p=",p_value,")")) +
  xlab("Mean PGA") + 
  ylab("mean Immune Clusters") + 
  xlim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/mICs_mPGA_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(summary(lm(Total_Events ~ Cont_Gleason, data = collected))$coefficients[2,4], digits = 3)

ggplot(collected, aes(x = Cont_Gleason, y = Total_Events)) + geom_point() + 
  ggtitle(paste0("Total Events v Gleason (p=",p_value,")")) +
  ylab("Total Events") + 
  xlab("Gleason") + 
  xlim(3,5) + 
  ylim(0,300) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_Total_Events_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(summary(lm(Subclonality ~ Cont_Gleason, data = collected))$coefficients[2,4], digits = 3)

ggplot(collected, aes(x = Cont_Gleason, y = Subclonality)) + geom_point() + 
  ggtitle(paste0("Subclonality v Gleason (p=",p_value,")")) +
  ylab("Subclonality") + 
  xlab("Gleason") + 
  xlim(3,5) + 
  ylim(0,1) + 
  geom_smooth(method = "lm") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_Subclonality_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(Cont_Gleason ~ TP53, colNArm(collected, sub = "TP53"))$p.value, digits = 3)
p_value_line_height = max(colNArm(collected, sub = "TP53")$Cont_Gleason, na.rm = T)*1.1

# Let's have a look at gleason and TP53
ggplot(colNArm(collected, sub = "TP53"), aes(x = TP53, y = Cont_Gleason)) + geom_boxplot() + 
  ggtitle(paste0("TP53 mutation status and Gleason")) +
  xlab("TP53 mutation?") + 
  ylab("Gleason") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height * 0.99, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height * 1.05, size = 5) +
  ylim(3,6) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_TP53_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(mPGA ~ TP53_clonal, colNArm(collected, sub = "TP53_clonal"))$p.value, digits = 3)
p_value_line_height = max(collected$mPGA)*1.1

ggplot(colNArm(collected, sub = "TP53_clonal"), aes(x = TP53_clonal, y = mPGA)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and mPGA")) +
  xlab("Clonal TP53?") + 
  ylab("Mean PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_clonal_PGA_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(max_PGA ~ TP53_clonal, colNArm(collected, sub = "TP53_clonal"))$p.value, digits = 3)
p_value_line_height = max(collected$max_PGA)*1.1

ggplot(colNArm(collected, sub = "TP53_clonal"), aes(x = TP53_clonal, y = max_PGA)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and max PGA")) +
  xlab("Clonal TP53?") + 
  ylab("Max PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_clonal_max_PGA_correlation.pdf", height = 6, width = 7)

ggplot(colNArm(collected, sub = "TP53_clonal"), aes(x = TP53_clonal, y = PsiT)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and PsiT")) +
  xlab("Clonal TP53?") + 
  ylab("Fitted Ploidy") + 
  ylim(1.5,4.5) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_clonal_PsiT_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(Total_Events ~ TP53_clonal, colNArm(collected, sub = "TP53_clonal"))$p.value, digits = 3)
p_value_line_height = max(colNArm(collected, sub = "TP53_clonal")$Total_Events)*1.1

# Let's have a look at events and TP53
ggplot(colNArm(collected, sub = "TP53_clonal"), aes(x = TP53_clonal, y = Total_Events)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and tree events")) +
  xlab("Clonal TP53?") + 
  ylab("Total Events") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height * 0.99, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height * 1.05, size = 5) +
  ylim(0,350) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_clonal_Events_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(Cont_Gleason ~ TP53_clonal, colNArm(collected, sub = "TP53_clonal"))$p.value, digits = 3)
p_value_line_height = max(colNArm(collected, sub = "TP53_clonal")$Cont_Gleason, na.rm = T)*1.1

# Let's have a look at gleason and TP53 clonality
ggplot(colNArm(collected, sub = "TP53_clonal"), aes(x = TP53_clonal, y = Cont_Gleason)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and Gleason")) +
  xlab("Clonal TP53?") + 
  ylab("Gleason") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height * 0.99, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height * 1.05, size = 5) +
  ylim(3,6) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Gleason_TP53_clonal_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(collected[collected$TP53_category=="Clonal","mPGA"], 
                        collected[collected$TP53_category=="Wild-type","mPGA"])$p.value, digits = 3)
p_value_line_height = max(colNArm(collected, sub = "TP53_category")$mPGA, na.rm = T)*1.1

ggplot(colNArm(collected, sub = "TP53_category"), aes(x = TP53_category, y = mPGA)) + geom_boxplot() + 
  ggtitle(paste0("Clonal TP53 mutation status and mPGA")) +
  xlab("Clonal TP53?") + 
  ylab("Mean PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 3, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 3, y = p_value_line_height * 0.99, xend = 3, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 2, y = p_value_line_height * 1.05, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/TP53_clonality_PGA_correlation.pdf", height = 6, width = 7)

# P Values
p_value = signif(t.test(collected[collected$DNA_damage_mutation=="Clonal","mPGA"], 
                        collected[collected$DNA_damage_mutation=="Wild-type","mPGA"])$p.value, digits = 3)
p_value_line_height = max(colNArm(collected, sub = "DNA_damage_mutation")$mPGA, na.rm = T)*1.1

# DNA damage mutation
ggplot(colNArm(collected, sub = "DNA_damage_mutation"), aes(x = DNA_damage_mutation, y = mPGA)) + geom_boxplot() + 
  ggtitle(paste0("DNA damage mutation (TP53/ATM/BRCA/CDK12/PALB2)\n and mPGA")) +
  xlab("") + 
  ylab("Mean PGA") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 3, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 3, y = p_value_line_height * 0.99, xend = 3, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 2, y = p_value_line_height * 1.05, size = 5) +
  ylim(0,1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/DNA_damage_clonality_PGA_correlation.pdf", height = 6, width = 7)

plot(collected[,c("mPGA", "Total_Events", "Number_subclonal_events")], pch = 20)

plots = lapply(c("mPGA", "Lossness", "Spearman", "Total_Events", "Number_subclonal_events", "Subclonality"), function(i) {
  len = length(unlist(strsplit(i, split = "_")))
  xl  = paste0(c(gsub("_", "\n", i), rep("\n", times = 3-len)), collapse = "")
  p = ggplot(collected, aes_string(x = 0, y = i)) + geom_violin() + geom_boxplot(width=0.1) + 
    xlab(xl) + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  return(p)
})

ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]], 
          labels = LETTERS[7:12],
          ncol = 6, nrow = 1)
ggsave("results/8_collect_metrics/Metric_distributions_genomics.png", height = 3, width = 10)



