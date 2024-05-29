# Analysis for associating chromosome arm changes with digital pathology

library(plyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(lmerTest)

# Functions
source("scripts/00_general_functions.R")

# Minimum PGA
min_pga = 0.01

# Get all files
call_files = list.files("results/0_min_het_best_fit_ploidy_search/calls/", 
                        pattern = "_cna_ploidy_search_calls.txt",
                        full.names = T)

# Comparing per sample image analysis data to per sample metrics
per_sample_metrics = read.table("results/2_metrics/metric_assessment_table.txt", stringsAsFactors = F, header = T)

# Remove normal samples and those without CNAs
per_sample_metrics = per_sample_metrics[!(per_sample_metrics$PGA < min_pga | per_sample_metrics$Assessment=="Failed" 
                                          | grepl("_N[0-9]{1,2}_", per_sample_metrics$Sample) | per_sample_metrics$Sample==""),]

# Read in gleason
gleason = read.csv("refs/GleasonSummaryClean.csv", 
                   stringsAsFactors = F, header = T)
gleason$Patient = unlist(lapply(strsplit(gleason$X, split = ""), function(i) paste0(i[1:6], collapse = "")))

# Read in the new per region data
gleason_regions = read.csv("refs/GleasonSummary_Regions.csv", 
                           stringsAsFactors = F, header = T)
gleason_regions$Patient = unlist(lapply(strsplit(gleason_regions$X, split = ""), function(i) paste0(i[1:6], collapse = "")))

# Get conversion reference table
look_up_table_slides = read.csv("refs/Forecast ID Reference.csv", 
                                stringsAsFactors = F, header = T)

# Filter for those we have sequencing data
sequenced_slides = look_up_table_slides[look_up_table_slides$Region!="",]
sequenced_slides$Region = gsub(" ", "", sequenced_slides$Region)

# Add things to use for identifying samples
sequenced_slides$Sample = paste0(sequenced_slides$Patient,"_",sequenced_slides$Region,"_DW1")
sequenced_slides$X = paste0(sequenced_slides$Patient,sequenced_slides$Slide)

# Extract only those we have sequenced slides
gleason = merge(gleason, sequenced_slides, by = "X")

# Change name format
gleason_regions$Sample = unlist(lapply(strsplit(gleason_regions$X, split = ""), 
                                       function(i) paste0(paste0(i[1:6], collapse = ""),"_",paste0(i[7:length(i)], collapse = ""),"_DW1")))

# Prepare for combination
gleason$Patient.x = NULL
colnames(gleason)[10] = "Patient"
gleason$Block = NULL
gleason$Slide = NULL
gleason$Region = NULL

# Gleason regions
gleason = rbind(gleason, gleason_regions)

# Add a continuous gleason score
gleason$Cont_Gleason = apply(gleason, 1, function(i) {
  
  sum(as.numeric(i[c("Gleason.3", "Gleason.4", "Gleason.5")])*c(3:5)) / sum(as.numeric(i[c("Gleason.3", "Gleason.4", "Gleason.5")]))
  
})

# Output the product
write.csv(gleason, "results/8_collect_metrics/Gleason_per_sample_digital_pathology.csv", quote = F, row.names = F)

# Subset only for samples with CNAs
call_files = call_files[gsub("_cna_ploidy_search_calls.txt", "", basename(call_files)) %in% per_sample_metrics$Sample]

# Read in the data
cn_data = lapply(call_files, read.table, header = T, stringsAsFactors = F)

# Read in some arm boundaries
arms = read.table("refs/chrArmBoundaries_hg38.txt", header = T)

chr_pos = extractChrRows(rownames(cn_data[[1]]))

chr_pos$chr = factor(chr_pos$chr, levels = c(1:22,"X","Y"))

# Function for assessing chromosome arm
getPQ = function(df, arm_pos) {
  
  psnqs = unlist(lapply(1:nrow(df), function(i) {
    
    r = df[i,]
    
    pos = (r$end - r$start) + r$start
    
    chr = r$chr
    
    b = arm_pos[which(arm_pos[,1] == chr),2]
    
    if(pos <  b) {a = "p"}
    if(pos >= b) {a = "q"}
    
    return(a)
    
  }))
  
  return(psnqs)
  
}

# Summarise the CNs of each arm
per_bin_arm = paste0("chr",chr_pos$chr,getPQ(chr_pos[,c(1:3)], arms))

# Make factor
per_bin_arm = factor(per_bin_arm, levels = unique(per_bin_arm), ordered = TRUE)

# Combine
cn_data = do.call(cbind, cn_data)

# Condense that per arm as median
median_cn_per_arm = apply(cn_data, 2, function(i) round(tapply(i, per_bin_arm, median)))

# Remove arms that aren't useful
median_cn_per_arm = median_cn_per_arm[!rownames(median_cn_per_arm) %in% c("chr21p", "chr22p"),]

# Get ploidy 
ploidy = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                    header = T, sep = "\t", stringsAsFactors = F)
colnames(ploidy)[1] = "Patient"

# Round it off
ploidy$Ploidy = round(ploidy$PsiT)

#################################################################################################################

# Get the gain/loss status
gain_cn_per_arm = do.call(cbind, lapply(1:ncol(median_cn_per_arm), function(col) {
  
  arm_cn  = median_cn_per_arm[,col]
  patient = unlist(strsplit(colnames(median_cn_per_arm)[col], split = "_"))[1]
  
  pld = ploidy[ploidy$Patient==patient,"Ploidy"]
  
  pld = rep(pld, nrow(median_cn_per_arm))
  
  pld[grepl("chrX|chrY", rownames(median_cn_per_arm))] = pld[grepl("chrX|chrY", rownames(median_cn_per_arm))]/2
  
  gains = arm_cn > pld
  
  return(gains)
  
}))

# Add column names
colnames(gain_cn_per_arm) = colnames(median_cn_per_arm)

# Run through the samples and gather genetic data
gain_cn_per_arm = as.data.frame(t(gain_cn_per_arm))
gain_cn_per_arm$Sample = rownames(gain_cn_per_arm)

#################################################################################################################

# Get the gain/loss status
loss_cn_per_arm = do.call(cbind, lapply(1:ncol(median_cn_per_arm), function(col) {
  
  arm_cn  = median_cn_per_arm[,col]
  patient = unlist(strsplit(colnames(median_cn_per_arm)[col], split = "_"))[1]
  
  pld = ploidy[ploidy$Patient==patient,"Ploidy"]
  
  pld = rep(pld, nrow(median_cn_per_arm))
  
  pld[grepl("chrX|chrY", rownames(median_cn_per_arm))] = pld[grepl("chrX|chrY", rownames(median_cn_per_arm))]/2
  
  losss = arm_cn < pld
  
  return(losss)
  
}))

# Add column names
colnames(loss_cn_per_arm) = colnames(median_cn_per_arm)

# Run through the samples and gather genetic data
loss_cn_per_arm = as.data.frame(t(loss_cn_per_arm))
loss_cn_per_arm$Sample = rownames(loss_cn_per_arm)

#################################################################################################################

# Format arms
arms = rownames(median_cn_per_arm)
arms = factor(arms, levels = unique(arms))

#################################################################################################################

loss_cn_per_arm$Sample = NULL
gain_cn_per_arm$Sample = NULL

gain_loss_cn_arms = gain_cn_per_arm

gain_loss_cn_arms = t(apply(gain_loss_cn_arms, 1, function(i) ifelse(i, "gain", "neutral")))

for(row in 1:nrow(loss_cn_per_arm)) {
  
  for(col in 1:ncol(loss_cn_per_arm)) {
    
    if(loss_cn_per_arm[row, col]) {gain_loss_cn_arms[row,col] = "loss"}
    
  }
  
}

gain_loss_cn_arms = as.data.frame(gain_loss_cn_arms)

for(col in 1:ncol(gain_loss_cn_arms)) {
  
  gain_loss_cn_arms[,col] = factor(gain_loss_cn_arms[,col], levels = c("neutral", "loss", "gain"))
  
}

gain_loss_cn_arms$Sample = rownames(gain_loss_cn_arms)

# Merge it together
gain_loss_cn_arms_gleason    = merge(gleason, gain_loss_cn_arms, by = "Sample")

# Calculate p values from mixed linear model
loss_ps = lapply(rownames(median_cn_per_arm), function(arm) {
  
  df_sub = gain_loss_cn_arms_gleason[gain_loss_cn_arms_gleason[,arm]!="gain",c("Patient", "Cont_Gleason", arm)]
  
  len = length(which(df_sub[,arm]=="loss"))
  
  if(len > 10) {
    
    lm = lmer(as.formula(paste0("Cont_Gleason ~ ",arm," + (1 | Patient)")), REML = TRUE, data = df_sub)

    p = summary(lm)$coefficients[2,5]
    m = summary(lm)$coefficient[2,"Estimate"]
    
    plt = ggplot(df_sub, aes_string(x = arm, y = "Cont_Gleason")) + geom_violin() + geom_boxplot(width=0.1) + 
      ggtitle(paste0("p=",signif(p, digits = 3)," m=",signif(m, digits = 3)))
    
    print(plt)
    
  } else {
    p = NA
    m = NA
  }
  
  out = list()
  
  out$p = p
  out$m = m
  out$n = len
  
  return(out)
  
})

names(loss_ps) = rownames(median_cn_per_arm)

# Calculate p values from mixed linear model
gain_ps = lapply(rownames(median_cn_per_arm), function(arm) {
  
  df_sub = gain_loss_cn_arms_gleason[gain_loss_cn_arms_gleason[,arm]!="loss",c("Patient", "Cont_Gleason", arm)]
  
  len = length(which(df_sub[,arm]=="gain"))
  
  if(len > 10) {
    
    lm = lmer(as.formula(paste0("Cont_Gleason ~ ",arm," + (1 | Patient)")), REML = TRUE, data = df_sub)
    
    p = summary(lm)$coefficients[2,5]
    m = summary(lm)$coefficient[2,"Estimate"]
    
    plt = ggplot(df_sub, aes_string(x = arm, y = "Cont_Gleason")) + geom_violin() + geom_boxplot(width=0.1) + 
      ggtitle(paste0("p=",signif(p, digits = 3)," m=",signif(m, digits = 3)))
    
    print(plt)
    
  } else {
    p = NA
    m = NA
  }
  
  out = list()
  
  out$p = p
  out$m = m
  out$n = len
  
  return(out)
  
})

names(gain_ps) = rownames(median_cn_per_arm)

# All of them together
res = data.frame(Type = rep(c("gain", "loss"), each = length(gain_ps)), 
                 Arm = rep(arms, times = 2), 
                 p = c(unlist(lapply(gain_ps, function(i) i$p)), unlist(lapply(loss_ps, function(i) i$p))),
                 m = c(unlist(lapply(gain_ps, function(i) i$m)), unlist(lapply(loss_ps, function(i) i$m))),
                 n = c(unlist(lapply(gain_ps, function(i) i$n)), unlist(lapply(loss_ps, function(i) i$n))))

# Calculate q value
res$p.adjust = p.adjust(p = res$p, method = "fdr")
res$pa_signif = res$p.adjust < 0.05

# Plot the results
plot_df = merge(res[res$Type=="gain",], res[res$Type=="loss",], by = "Arm")
p = ggplot(plot_df, aes(x = Arm, y = -1*log(p.adjust.x, base = 10))) + geom_bar(aes(fill = "gains"), stat = "identity") + 
  geom_bar(aes(y=log(p.adjust.y, base = 10), fill="losses"), stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") + 
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  ggtitle("Association of Chromosome Arm Changes with Gleason Difference") +
  scale_y_continuous(name = "log10(p.adjust)", limits = c(-6,6), breaks = seq(-6,6,by=2), labels = -abs(seq(-6,6,by=2))) + 
  xlab("") + geom_hline(yintercept = c(log(0.05, base = 10), -1*log(0.05, base = 10)), lty = "dotted") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p, filename = "results/8_collect_metrics/Association_chromosome_arms_gleason_difference.png", height = 5, width = 8)

# Volcano plot
p = ggplot(res, aes(x = m, y = -1 * log(p.adjust, base = 10), col = Type, 
                    label = paste0(ifelse(Type=="gain", "+", "-"),gsub("chr", "", Arm)))) + 
  geom_text() +
  geom_hline(yintercept = -1*log(0.05, base = 10), lty = "dotted") + theme(legend.position = "none") +
  ylab(expression(-log[10]*italic("p"))) +
  xlab("Continuous Gleason change") + 
  ggtitle("Per arm change in Gleason") +
  scale_color_manual(values=c("loss"="#3c73a8","gain"="#ec2d01")) + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave(p, filename = "results/8_collect_metrics/Volcano_plot_arms_gleason_difference.png", height = 3, width = 4)

# Read in the arm data
tsgog_score = read.csv("refs/Davoli_2013_Cell_TSG-OG_scores.csv", stringsAsFactors = F)
tsgog_score$Charm_TSG_OG_score = tsgog_score$Charm_TSG_OG_score * -1
tsgog_score$Arm                = paste0("chr", tsgog_score$Arm)

# Add the scores
res = merge(res, tsgog_score, by = "Arm")

# Add directionality too
res$pa_signif_direction = ifelse(res$pa_signif, yes = paste0(res$pa_signif,"_",res$m>0), res$pa_signif)

# Remove NAs
res_omit = na.omit(res)

# Let's do some t.tests
res_gain = res_omit[res_omit$Type == "gain" & res_omit$pa_signif_direction %in% c("FALSE", "TRUE_TRUE"),]
res_loss = res_omit[res_omit$Type == "loss" & res_omit$pa_signif_direction %in% c("FALSE", "TRUE_TRUE"),]

# Positive v none pvalues
p1 = signif(t.test(Charm_TSG_OG_score ~ pa_signif_direction, res_gain, alternative = "less")$p.value, digits = 3)
p2 = signif(t.test(Charm_TSG_OG_score ~ pa_signif_direction, res_loss, alternative = "greater")$p.value, digits = 3)

# Make as dataframe
ps      = data.frame(pa_signif_direction = c("TRUE_FALSE"), 
                     Charm_TSG_OG_score = c(6,4.5), Type = c("gain", "loss"))
sg      = data.frame(x = c(1,1,3,1,1,3), y = c(5.5,5.5,5.5,4,4,4), 
                     xend = c(3,1,3,3,1,3), yend = c(5.5,5.3,5.3,4,3.8,3.8),
                     Type = c(rep("gain", times = 3), rep("loss", times = 3)))

# Associate to TSG-OG score
p = ggplot(res_omit, aes(x = pa_signif_direction, y = Charm_TSG_OG_score)) + 
  geom_violin(width=0.5) + geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_jitter(width = 0.2) + 
  ylim(-5,6) +
  facet_grid(cols = vars(Type)) + xlab("Significant Association with Gleason") + ylab("TSG-OG score") + 
  geom_text(data = ps, label = paste0("p=",c(p1, p2))) +
  geom_segment(data = sg, aes(x = x, xend = xend, y = y, yend = yend)) +
  scale_x_discrete(labels=c("None", "Negative", "Positive"))
ggsave(p, filename = "results/8_collect_metrics/Chromosome_Arm_Change_TSG_OG_scores.png", height = 5, width = 5)

#################################################################################################################

# Read in Dropbox
morisita = read.csv("refs/Morisita_Cancer_Regions.csv", 
                    stringsAsFactors = F)

# Edit patient/sample names
morisita$Sample = convertPathSampleName(morisita$SampleID)

# Add patient id
morisita$Patient = unlist(lapply(strsplit(morisita$Sample, split = ""), function(i) paste0(i[1:6], collapse = "")))

# Subset for only passing samples
gain_loss_cn_arms = gain_loss_cn_arms[gain_loss_cn_arms$Sample %in% per_sample_metrics$Sample,]

# Add it on
gain_loss_cn_arms_timorisita = merge(morisita, gain_loss_cn_arms, by = "Sample")

# Mixed effects model for gleason morisita 
loss_ps = lapply(rownames(median_cn_per_arm), function(arm) {
  
  df_sub = gain_loss_cn_arms_timorisita[gain_loss_cn_arms_timorisita[,arm]!="gain",c("Patient", "Morisita", arm)]
  
  len = length(which(df_sub[,arm]=="loss"))
  
  if(len > 10) {
    
    lm = lmer(as.formula(paste0("Morisita ~ ",arm," + (1 | Patient)")), REML = TRUE, data = df_sub)
    
    p = summary(lm)$coefficients[2,5]
    m = summary(lm)$coefficient[2,"Estimate"]
    
  } else {
    p = NA
    m = NA
  }
  
  out = list()
  
  out$p = p
  out$m = m
  out$n = len
  
  return(out)
  
})

# Format output
names(loss_ps) = rownames(median_cn_per_arm)

# Mixed effects model for Gleason Morisita
gain_ps = lapply(rownames(median_cn_per_arm), function(arm) {
  
  df_sub = gain_loss_cn_arms_timorisita[gain_loss_cn_arms_timorisita[,arm]!="loss",c("Patient", "Morisita", arm)]
  
  len = length(which(df_sub[,arm]=="gain"))
  
  if(len > 10) {
    
    lm = lmer(as.formula(paste0("Morisita ~ ",arm," + (1 | Patient)")), REML = TRUE, data = df_sub)
    
    p = summary(lm)$coefficients[2,5]
    m = summary(lm)$coefficient[2,"Estimate"]
    
  } else {
    p = NA
    m = NA
  }
  
  out = list()
  
  out$p = p
  out$m = m
  out$n = len
  
  return(out)
  
})

# Format output
names(gain_ps) = rownames(median_cn_per_arm)

# All of them together
res = data.frame(Type = rep(c("gain", "loss"), each = length(gain_ps)), 
                 Arm = rep(arms, times = 2), 
                 p = c(unlist(lapply(gain_ps, function(i) i$p)), unlist(lapply(loss_ps, function(i) i$p))),
                 m = c(unlist(lapply(gain_ps, function(i) i$m)), unlist(lapply(loss_ps, function(i) i$m))),
                 n = c(unlist(lapply(gain_ps, function(i) i$n)), unlist(lapply(loss_ps, function(i) i$n))))

# Calculate q value
res$p.adjust = p.adjust(p = res$p, method = "fdr")
res$pa_signif = res$p.adjust < 0.05

p = ggplot(res, aes(x = m, y = -1 * log(p.adjust, base = 10), col = Type, 
                    label = paste0(ifelse(Type=="gain", "+", "-"),gsub("chr", "", Arm)))) + geom_text() +
  geom_hline(yintercept = -1*log(0.05, base = 10), lty = "dotted") + theme(legend.position = "none") +
  ylab(expression(-log[10]*italic("p"))) +
  xlab("Tumour-Immune Morisita change") + 
  scale_color_manual(values=c("loss"="#3c73a8","gain"="#ec2d01")) + 
  ggtitle("Per arm change in TI Morisita") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(p, filename = "results/8_collect_metrics/Volcano_plot_arms_tumour_immune_morisita_difference.png", height = 3, width = 4)

# Show the chr6p result
arm = "chr6p"

df_sub = gain_loss_cn_arms_timorisita[gain_loss_cn_arms_timorisita[,arm]!="gain",c("Patient", "Morisita", arm)]

p_value_line_height = 0.98
p_value = signif(res[res$Type=="loss" & res$Arm=="chr6p","p"], digits = 2)

plt = ggplot(df_sub, aes_string(x = arm, y = "Morisita")) + ylab("Tumour-Immune Morisita") + geom_violin() + 
  geom_segment(aes(x = 1, y = p_value_line_height, xend = 2, yend = p_value_line_height)) +
  geom_segment(aes(x = 1, y = p_value_line_height * 0.99, xend = 1, yend = p_value_line_height)) +
  geom_segment(aes(x = 2, y = p_value_line_height * 0.99, xend = 2, yend = p_value_line_height)) +
  annotate("text", label = paste0("p=",p_value), x = 1.5, y = p_value_line_height * 1.1, size = 5) +
  scale_y_continuous(breaks=c(0,0.5,1), limits = c(0,1.2)) +
  geom_boxplot(width=0.1) + scale_x_discrete(labels=c("Baseline", "Loss")) + xlab("Chromosome 6p")

ggsave(plt, filename = "results/8_collect_metrics/TI_Morisita_chr6_association.png", height = 3, width = 4)

#################################################################################################################

# Show immune infiltration per chromosome 6p status
cell_class = read.csv("CellClassificationSummary_Regions.csv",
                      stringsAsFactors = F)

# Convert sample naming
cell_class$Sample = convertPathSampleName(cell_class$X)

# Merge
calls_cell_class = merge(gain_loss_cn_arms_timorisita, cell_class, by = "Sample")

calls_cell_class$chr6p_loss_patient = calls_cell_class$Patient %in% unique(calls_cell_class$Patient[calls_cell_class$chr6p == "loss"])

lm = lmer(as.formula(paste0("Immune.. ~ chr6p + (1 | Patient)")), REML = TRUE, data = calls_cell_class[calls_cell_class$chr6p!="gain",])

p = signif(summary(lm)$coefficients[2,5], digits = 3)
m = summary(lm)$coefficient[2,"Estimate"]

plt = ggplot(calls_cell_class, aes(x = chr6p, y = Immune..)) + geom_boxplot() + 
  ylab("Immune %") + ylim(0,75) +
  geom_segment(aes(x = 1, y = 60, xend = 2, yend = 60)) +
  geom_segment(aes(x = 1, y = 60 - 2, xend = 1, yend = 60)) +
  geom_segment(aes(x = 2, y = 60 - 2, xend = 2, yend = 60)) +
  annotate("text", label = paste0("p=",p), x = 1.5, y = 60 + 5, 
           size = 3) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(plt, filename = "results/8_collect_metrics/Immune_infiltration_chr6_loss.png", height = 3, width = 5)
