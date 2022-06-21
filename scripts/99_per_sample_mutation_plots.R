# Per sample mutations plots

library(plyr)
library(ggplot2)
library(cowplot)
library(lmerTest)
library(reshape2)

# Functions
source("scripts/00_general_functions.R")

# Minimum PGA
min_pga = 0.01

# Targeted data repo
targeted_repo = "target_repo/"

# Comparing per sample image analysis data to per sample metrics
per_sample_metrics = read.table("results/2_metrics/metric_assessment_table.txt", stringsAsFactors = F, header = T)

# Remove low CNA samples and normals
per_sample_metrics = per_sample_metrics[!(per_sample_metrics$PGA < min_pga | per_sample_metrics$Assessment=="Failed" 
                                          | grepl("_N[0-9]{1,2}_", per_sample_metrics$Sample) | per_sample_metrics$Sample==""),]

# Read in the mutation data
per_sample_mutation = read.table(paste0(targeted_repo,"results/mutation_calling/Per_sample_gene_mutation_status.txt",)
                                 sep = "\t", header = T)

# Reformat the names
per_sample_mutation$Sample = paste0(per_sample_mutation$Patient,"_",per_sample_mutation$Sample)
per_sample_mutation$Patient = NULL

# Merge it
per_sample_data = merge(per_sample_metrics, per_sample_mutation, by = "Sample")

# We don't need NAs
per_sample_data = per_sample_data[!is.na(per_sample_data$AKT1),]

# Do mixed modelling
res2b = lmer(PGA ~ TP53 + (1 | Patient), 
             REML = TRUE, data = per_sample_data)
summary(res2b)

# P Values
p_value = signif(summary(res2b)$coefficients[2,5], digits = 2)
p_value_line_height = max(per_sample_data$PGA)*1.1

# Plot it for show
ggplot(per_sample_data, aes(x=TP53, y=PGA)) + ylim(0,1) + geom_violin() + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  annotate("text", label = paste0("n=",table(per_sample_data$TP53)), x = 1:2, y = c(0.9,0.9), 
           size = 5) +
  xlab("TP53 mutated") + geom_boxplot(width=0.1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Sample_based_MLM_TP53_mutation.png", height = 3, width = 3.5)

# DNA damage
per_sample_data$DNA_damage = per_sample_data$CDK12 | per_sample_data$ATM | per_sample_data$TP53 | per_sample_data$BRCA1 | per_sample_data$BRCA2 | per_sample_data$PALB2

# Do mixed modelling
res2b = lmer(PGA ~ DNA_damage + (1 | Patient), 
             REML = TRUE, data = per_sample_data)
summary(res2b)

# P Values
p_value = signif(summary(res2b)$coefficients[2,5], digits = 2)
p_value_line_height = max(per_sample_data$PGA)*1.1

ggplot(per_sample_data, aes(x=DNA_damage, y=PGA)) + ylim(0,1) + geom_violin() + 
  ggtitle("PGA and DNA damage mutation, mixed linear model\n(CDK12, ATM, TP53, BRCA1/2, PALB2)") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.01, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.01, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.075, size = 5) +
  annotate("text", label = paste0("n=",table(per_sample_data$DNA_damage)), x = 1:2, y = c(0.9,0.9), 
           size = 5) +
  xlab("DNA damage mutated") + geom_boxplot(width=0.1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Sample_based_MLM_DNA_damage_mutation.pdf", height = 6, width = 7)

#################################
##### Add Gleason stuff too #####
#################################

# Read in the gleason per sample data
gleason = read.csv("results/8_collect_metrics/Gleason_per_sample_digital_pathology.csv", 
                   header = T, stringsAsFactors = F)

gleason$Patient = NULL

# Merge retaining left side
per_sample_data = merge(per_sample_data, gleason, by = "Sample", all.x = T)

# Do mixed modelling
res2b = lmer(Cont_Gleason ~ TP53 + (1 | Patient), 
             REML = TRUE, data = per_sample_data)
summary(res2b)

# P Values
p_value = signif(summary(res2b)$coefficients[2,5], digits = 2)
p_value_line_height = max(per_sample_data$Cont_Gleason, na.rm = T)*1.1

# Compare TP53 mutation status to continuous Gleason
ggplot(per_sample_data, aes(x=TP53, y=Cont_Gleason)) + ylim(3,5.8) + geom_violin() + 
  ylab("Continuous Gleason") +
  ggtitle("Continuous Gleason and TP53 mutation, mixed linear model") + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.05, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.05, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height + 0.2, size = 5) +
  annotate("text", label = paste0("n=",table(na.omit(per_sample_data[,c("TP53", "Cont_Gleason")])$TP53)), x = 1:2, y = c(5.2,5.2), 
           size = 5) +
  xlab("TP53 mutated") + geom_boxplot(width=0.1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Sample_based_MLM_Cont_Gleason_fraction_TP53.pdf", height = 6, width = 7)

ggplot(per_sample_data, aes(x=TP53, y=Cont_Gleason)) + ylim(3,5.8) + geom_violin() +
  ylab("Continuous Gleason") +
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height - 0.05, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height - 0.05, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",format(p_value, scientific = FALSE)), x = 1.5, y = p_value_line_height + 0.2, size = 5) +
  annotate("text", label = paste0("n=",table(na.omit(per_sample_data[,c("TP53", "Cont_Gleason")])$TP53)), x = 1:2, y = c(5.2,5.2),
  size = 5) +
  xlab("TP53 mutated") + geom_boxplot(width=0.1) + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/8_collect_metrics/Sample_based_MLM_Cont_Gleason_fraction_TP53.png", height = 3, width = 4)
