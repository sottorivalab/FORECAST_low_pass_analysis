# Compare per sample digital pathology results to genomics

library(plyr)
library(ggplot2)
library(cowplot)

# Functions
source("scripts/00_general_functions.R")

# Minimum PGA
min_pga = 0.01

# Comparing per sample image analysis data to per sample metrics
per_sample_metrics = read.table("results/2_metrics/metric_assessment_table.txt", stringsAsFactors = F, header = T)

# Remove low CNA samples
per_sample_metrics = per_sample_metrics[!(per_sample_metrics$PGA < min_pga | per_sample_metrics$Assessment=="Failed" 
                                          | grepl("_N[0-9]{1,2}_", per_sample_metrics$Sample) | per_sample_metrics$Sample=="IM1017_B2_DW1"),]

# Read in tumour morisita
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
collected = merge(per_sample_metrics, morisita, by = "Sample")
collected = merge(collected, lymphocyte, by = "Sample")
collected = merge(collected, cell_class, by = "Sample")

# Rounded ploidy
collected$Round_ploidy = as.factor(round(collected$Ploidy/2)*2)

# Cut
collected$PGA_cut = cut(collected$PGA, breaks = c(0,quantile(collected$PGA, probs = 1:3/4),1))

pv1 = signif(wilcox.test(collected$Morisita[collected$PGA_cut==levels(collected$PGA_cut)[3]], 
                         collected$Morisita[collected$PGA_cut %in% c(levels(collected$PGA_cut)[1], 
                                                                     levels(collected$PGA_cut)[2])])$p.value, 
             digits = 2)
pv2 = format(wilcox.test(collected$Morisita[collected$PGA_cut==levels(collected$PGA_cut)[4]], 
                         collected$Morisita[collected$PGA_cut %in% c(levels(collected$PGA_cut)[1], 
                                                                     levels(collected$PGA_cut)[2])])$p.value, 
             digits = 2, 
             scientific = F)

# Cutting PGA up into ranges
ggplot(collected, aes(x = PGA_cut, y = Morisita)) + geom_boxplot() + 
  xlab("PGA Quartiles") + ylab("Tumour-Immune Morisita") + 
  scale_x_discrete(labels=c("1st", "2nd", "3rd", "4th")) +
  geom_segment(aes(x = 1.5, y = 0.98, xend = 3, yend = 0.98)) +
  geom_segment(aes(x = 1.5, y = 0.98 * 0.99, xend = 1.5, yend = 0.98)) +
  geom_segment(aes(x = 3, y = 0.98 * 0.99, xend = 3, yend = 0.98)) +
  geom_segment(aes(x = 1.5, y = 1.1, xend = 4, yend = 1.1)) +
  geom_segment(aes(x = 1.5, y = 1.1 * 0.99, xend = 1.5, yend = 1.1)) +
  geom_segment(aes(x = 4, y = 1.1 * 0.99, xend = 4, yend = 1.1)) +
  annotate("text", label = paste0("p=",pv1), x = 2.25, y = 0.98 * 1.06, size = 4) +
  annotate("text", label = paste0("p=",pv2), x = 2.75, y = 1.1 * 1.06, size = 4) +
  scale_y_continuous(breaks=c(0,0.5,1), limits = c(0,1.2)) +
  theme_cowplot()
ggsave("results/7_collect_metrics/TI_Morisita_PGA_quartiles.png", height = 3, width = 4)
