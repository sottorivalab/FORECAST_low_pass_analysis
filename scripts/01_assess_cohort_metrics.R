# Plot metrics of the cohort
library(ggplot2)
library(cowplot)

# Functions
source("scripts/00_general_functions.R")

# Where is the low pass data?
pipeline_root = "~/remotes/forecast_low_pass/"

# Read in the manual assessment of the lpWGS
low_pass_assessment = read.table("refs/manual_assessment_categories.txt",
                                 header = T)

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Pipeline metrics files
per_patient_samples = lapply(patients, FORECAST_list.files, pattern = "_metrics.txt")

# Read in the data
per_patient_metrics = lapply(per_patient_samples, FORECAST_read.files)

# Make per_patient matrix
per_patient_met_mat = lapply(per_patient_metrics, function(p) do.call(rbind, p))

# Make a cohort level matrix
per_patient_met_mat_pipe = do.call(rbind, per_patient_met_mat)

# Assess quality
per_patient_met_mat_pipe$Quality = "Good"
per_patient_met_mat_pipe$Quality[per_patient_met_mat_pipe$Within_segment_variance > 0.09] = "Failed"
per_patient_met_mat_pipe$Quality[per_patient_met_mat_pipe$Within_segment_variance > 0.06 &
                                 per_patient_met_mat_pipe$Within_segment_variance <=0.09] = "Bad"

# Plot number of reads and variance
distribute_annotation_x = sample(c(0.33,0.5,2,3), 
                                 size = count(per_patient_met_mat_pipe$Quality!="Good"), 
                                 replace = T)
distribute_annotation_y = distribute_annotation_x

ggplot(per_patient_met_mat_pipe, aes(x = Read_Count, y = Within_segment_variance, col = Quality)) + 
  geom_point() + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + 
  xlab("Binned Reads") + ylab("Variance within segments") + theme_cowplot() + 
  annotate("segment",
           x = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Read_Count"],
           xend = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Read_Count"] * 
             distribute_annotation_x,
           y = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Within_segment_variance"], 
           yend = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Within_segment_variance"] * 
             distribute_annotation_y,
           colour = "black", size=0.5, alpha=1) +
  annotate("text", 
           size = 3,
           x = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Read_Count"] * 
             (distribute_annotation_x * ifelse(distribute_annotation_x < 1, yes = 0.8, no = 1.2)), 
           y = per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Within_segment_variance"] * 
             (distribute_annotation_y* ifelse(distribute_annotation_y < 1, yes = 0.8, no = 1.2)),
           label = gsub("DW1", "", 
                        gsub("_", " ", 
                             paste0(per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Patient"],
                                    "_",per_patient_met_mat_pipe[per_patient_met_mat_pipe$Quality!="Good","Region"]))))
ggsave("results/2_metrics/Segment_variance_coverage_QC.pdf", height = 6, width = 7)

# Need a sample column
per_patient_met_mat_pipe$Sample = paste0(per_patient_met_mat_pipe$Patient,"_",per_patient_met_mat_pipe$Region)

# Combine with manual assessment
per_patient_met_mat_pipe = merge(per_patient_met_mat_pipe, low_pass_assessment, by = "Sample")

# Compare to manual assessment
ggplot(per_patient_met_mat_pipe, aes(x = Read_Count, y = Within_segment_variance, col = Assessment)) + 
  scale_colour_manual(values = c("#15b01a", "#fcb001", "#e50000", "#0343df")) +
  geom_point() + scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') + 
  xlab("Binned Reads") + ylab("Variance within segments") + theme_cowplot()
ggsave("results/2_metrics/Segment_variance_coverage_QC_manual_assessment.pdf", height = 6, width = 7)


# Get all files
metr_files = list.files("results/0_min_het_best_fit_ploidy_search/calls/", 
                        pattern = "_cna_ploidy_search_metrics.txt",
                        full.names = T)

# Find the metrics files
per_patient_samples = lapply(patients, function(p) metr_files[grep(p, x = metr_files)])

# Read in the data
per_patient_metrics = lapply(per_patient_samples, FORECAST_read.files)

# Make per_patient matrix
per_patient_met_mat = lapply(per_patient_metrics, function(p) do.call(rbind, p))

# Make a cohort level matrix
per_patient_met_mat = do.call(rbind, per_patient_met_mat)

# Add patient attribute
per_patient_met_mat$Patient = extractPatientName(per_patient_met_mat$Sample)

# Combine with manual assessment
per_patient_met_mat = merge(per_patient_met_mat, low_pass_assessment, by = "Sample")

# What is the mean after removing the zeros?
noZeroMeanPGA = mean(removeZero(per_patient_met_mat$PGA))

# Make a ggplot
ggplot(data = per_patient_met_mat, aes(x = PGA)) + 
  geom_histogram(bins = round(diff(range(per_patient_met_mat$PGA)) * 100), color = "black", fill = "gray") + 
  geom_vline(aes(xintercept = noZeroMeanPGA), linetype = "dashed", size = 0.6) + 
  ggtitle(paste0("PGA of all samples in the cohort (mean PGA = ",round(noZeroMeanPGA, digits = 2),")")) +
  ylab("Count") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/2_metrics/Histogram_all_PGA_cohort.pdf", height = 6, width = 7)

# Calculate mean PGA
mean_PGA = do.call(rbind, lapply(patients, function(pat) {
  
  mPGA = mean(removeZero(per_patient_met_mat[per_patient_met_mat$Patient==pat,"PGA"]))
  
  data.frame(Patient = pat, Mean_PGA = mPGA)
  
}))

# Create mean
mmPGA = mean(mean_PGA$Mean_PGA, na.rm = T)

# Add zeros in
mean_PGA$Mean_PGA[is.na(mean_PGA$Mean_PGA)] = 0

# Make a ggplot
ggplot(data = mean_PGA, aes(x = Mean_PGA)) + 
  geom_histogram(bins = round(diff(range(mean_PGA$Mean_PGA, na.rm = T)) * 100), color = "black", fill = "gray") + 
  geom_vline(aes(xintercept = mmPGA), linetype = "dashed", size = 0.6) +
  geom_vline(aes(xintercept = quantile(Mean_PGA, probs = 2/3)), linetype = "dotted", size = 0.6) +
  ggtitle(paste0("Mean PGA of each patient in the cohort (mean PGA = ",round(mmPGA, digits = 2),")")) +
  ylab("Count") +
  xlab("Mean PGA") + 
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/2_metrics/Histogram_mean_PGA_cohort.pdf", height = 6, width = 7)

# See if purity and PGA correlate
ggplot(data = per_patient_met_mat, aes(x = Purity, y = PGA)) + 
  geom_point() +
  ggtitle(paste0("Correlation of PGA and purity")) +
  ylab("PGA") +
  xlab("Purity") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/2_metrics/Purity_PGA_correlation.pdf", height = 6, width = 7)

# Write out results including the assessment of pass and fail
write.table(per_patient_met_mat, file = "results/2_metrics/metric_assessment_table.txt", sep = "\t", quote = F, 
            row.names = F)

# Get ploidy 
ploidy = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                    header = T, sep = "\t", stringsAsFactors = F)
colnames(ploidy)[1] = "Patient"

# Round it off
ploidy$Ploidy = round(ploidy$PsiT)

# Compare to mean PGA
mean_PGA_ploidy = merge(mean_PGA, ploidy, by = "Patient")

# Compare mean PGA to ploidy
ggplot(mean_PGA_ploidy, aes(x = as.factor(Ploidy), y = Mean_PGA)) + ylim(0,1) + 
  ggtitle(paste0("Mean PGA of Patient and inferred Ploidy")) +
  xlab("Ploidy") + 
  ylab("Mean PGA") + 
  geom_boxplot() + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/2_metrics/Ploidy_PGA_correlation.pdf", height = 6, width = 7)





