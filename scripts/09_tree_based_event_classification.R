# A script for classifying events by phylogenetic history and making summary plots

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(GenomicRanges)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38

library(phytools)
library(usedist)
library(adephylo)
library(ComplexHeatmap)
library(circlize)
options(scipen = 999)

# Functions
source("scripts/00_general_functions.R")

# Roots of data
medicc_root   = "~/remotes/medicc_trees/"
pipeline_root = "~/remotes/forecast_low_pass/"

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Read in the patient ploidies
patient_ploidies = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                              header = T, stringsAsFactors = F)
patient_ploidies$Ploidy = round(patient_ploidies$PsiT)

# Loop them
trees = lapply(patients, function(patient) {
  
  tree = try(read.tree(paste0(medicc_root,"forecast/",patient,"_regions_final_tree.new")))
  
  return(tree)
  
})

# Loop them
cnas = lapply(patients, function(patient) {
  
  res = try(read.table(paste0(medicc_root,"forecast/",patient,"_regions_final_cn_profiles.tsv"), 
                        header = T, stringsAsFactors = F))
  
  return(res)
  
})

# Tree based cna changes
cnas_tree_based = lapply(1:length(patients), function(i) {
  
  tree = trees[[i]]
  cna  = cnas[[i]]
  pat  = patients[[i]]
  
  if(class(tree)=="phylo") {
  
    breakdown = data.frame(A = tree$edge[,1],
                           B = tree$edge[,2],
                           Edge_Length = tree$edge.length, 
                           A_Label = c(tree$tip.label, tree$node.label)[tree$edge[,1]], 
                           B_Label = c(tree$tip.label, tree$node.label)[tree$edge[,2]],
                           stringsAsFactors = F)
    
    breakdown$A_Label[breakdown$A_Label==""] = "diploid"
    
    breakdown = breakdown[!(breakdown$A_Label == "diploid" & breakdown$B_Label == "diploid"),]
    
    clonal_branch = breakdown[breakdown$A_Label=="diploid",]
    
    root = cna[cna$sample_id==clonal_branch$A_Label,]
    
    mrca = cna[cna$sample_id==clonal_branch$B_Label,]
    
    clonal_cnas = mrca$cn - (root$cn * patient_ploidies[patient_ploidies$Case==pat,"Ploidy"] / 2)
    
    intermediate_branch = breakdown[grepl("internal_", breakdown$B_Label) & breakdown$A_Label!="diploid",]
    
    if(nrow(intermediate_branch) > 0) {
    
      intermediate_cnas = apply(intermediate_branch, 1, function(i) cna[cna$sample_id==i["B_Label"],"cn"] - cna[cna$sample_id==i["A_Label"],"cn"])
      
    } else {intermediate_cnas = cbind(rep(0, times = nrow(root)))}
    
    tip_branch = breakdown[!grepl("internal_", breakdown$B_Label),]
    
    tip_cnas = apply(tip_branch, 1, function(i) cna[cna$sample_id==i["B_Label"],"cn"] - cna[cna$sample_id==i["A_Label"],"cn"])
    
    cna_back = cna[1:nrow(root),2:4]
    
    res = data.frame(cna_back,
                     clonal_loss = apply(cbind(clonal_cnas), 1, function(i) as.integer(any(i < 0))),
                     clonal_gain = apply(cbind(clonal_cnas), 1, function(i) as.integer(any(i > 0))),
                     intermediate_loss = apply(intermediate_cnas, 1, function(i) as.integer(any(i < 0))),
                     intermediate_gain = apply(intermediate_cnas, 1, function(i) as.integer(any(i > 0))),
                     tip_loss = apply(tip_cnas, 1, function(i) as.integer(any(i < 0))),
                     tip_gain = apply(tip_cnas, 1, function(i) as.integer(any(i > 0))))
    
  } else {res = NULL}
  
  return(res)

})

# Get all files
call_files = list.files("results/0_min_het_best_fit_ploidy_search/calls/", 
                        pattern = "_cna_ploidy_search_calls.txt",
                        full.names = T)

# Read in the data
per_patient_calls = lapply(call_files[1], FORECAST_read.files)

# Make windows object from first sample
windows       = extractChrRows(rownames(per_patient_calls[[1]][[1]]))
windows$name  = paste0("win",1:nrow(windows))
windows_gr     = GRanges(seqnames = paste0("chr", windows$chr), ranges = IRanges(start = windows$start, end = windows$end))

# Expand out per patient
cnas_expanded = lapply(cnas_tree_based, function(mat) {
  
  if(!is.null(mat)) {
  
    mat_gr = GRanges(seqnames = mat$chrom, ranges = IRanges(start = mat$start, end = mat$end))
    
    res    = mat[subjectHits(findOverlaps(windows_gr, mat_gr)),4:ncol(mat)]
    
    rownames(res) = NULL
    
  } else {res = NULL}
  
  return(res)
  
})

# Collect the summaries
cnas_expanded = cnas_expanded[unlist(lapply(cnas_expanded, function(i) !is.null(i)))]

# Summarise for cohort
forecast_tree_cnas = Reduce('+', cnas_expanded)

# Collect
data_per_window_summary = cbind(windows, forecast_tree_cnas)

# What's the middle?
chr_centers = 
  filter(D3GB::GRCh38.bands) %>%
  split(., .$chr) %>% 
  sapply(function(x) max(x$end)/2)

centromer_position = 
  filter(D3GB::GRCh38.bands, score=="acen") %>%
  split(., .$chr) %>% 
  sapply(function(x) mean(range(c(x$end,x$start))))

ylm = 120

plot_cna_recurence=
  data_per_window_summary %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)]) %>%
  mutate(chr=factor(chr, c(1:22,"X","Y"), ordered=TRUE)) %>%
  ggplot(aes(x=pos, y=clonal_gain+intermediate_gain+tip_gain)) + 
  geom_area(aes(alpha="tip", fill="gains")) +
  geom_area(aes(y=clonal_gain+intermediate_gain, alpha="intermediate", fill="gains")) +
  geom_area(aes(y=clonal_gain, alpha="clonal", fill="gains")) + 
  geom_area(aes(y=-clonal_loss-intermediate_loss-tip_loss, alpha="tip", fill="losses")) + 
  geom_area(aes(y=-clonal_loss-intermediate_loss, alpha="intermediate", fill="losses")) + 
  geom_area(aes(y=-clonal_loss, alpha="clonal", fill="losses")) + 
  geom_hline(yintercept=0) + 
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers), 
             aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(.~chr, scales="free_x", space="free_x") + 
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) + 
  scale_alpha_manual(values=c("clonal"=1, "intermediate"=0.5, "tip"=0.2)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  scale_x_continuous(breaks=c(0), labels="") + 
  scale_y_continuous(breaks=seq(-ylm, ylm, by=10), 
                     labels=abs(seq(-ylm, ylm, by=10)), 
                     limits=c(-ylm,ylm)) + 
  xlab("") + ylab("Number of cases") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}

# Plot it!
plot(plot_cna_recurence)

# Save the output
ggsave("results/4_clonality_frequency/genome_wide_clonality_tree_based_frequency.png", plot_cna_recurence, width=11, height=4.2)

ylm = 60

plot_cna_recurence=
  data_per_window_summary %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)]) %>%
  mutate(chr=factor(chr, c(1:22,"X","Y"), ordered=TRUE)) %>%
  ggplot(aes(x=pos, y=clonal_gain)) + 
  geom_area(aes(alpha="clonal", fill="gains")) +
  geom_area(aes(y=-clonal_loss, alpha="clonal", fill="losses")) + 
  geom_hline(yintercept=0) + 
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers), 
             aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(.~chr, scales="free_x", space="free_x") + 
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) + 
  scale_alpha_manual(values=c("clonal"=1, "intermediate"=0.5, "tip"=0.2)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  scale_x_continuous(breaks=c(0), labels="") + 
  scale_y_continuous(breaks=seq(-ylm, ylm, by=20), 
                     labels=abs(seq(-ylm, ylm, by=20)), 
                     limits=c(-ylm,ylm)) + 
  xlab("") + ylab("Number of cases") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}

# Save the output
ggsave("results/4_clonality_frequency/genome_wide_clonality_clonal_tree_based_frequency.png", plot_cna_recurence, width=11, height=2.1)

ylm = 60

plot_cna_recurence=
  data_per_window_summary %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)]) %>%
  mutate(chr=factor(chr, c(1:22,"X","Y"), ordered=TRUE)) %>%
  ggplot(aes(x=pos, y=intermediate_gain)) + 
  geom_area(aes(alpha="intermediate", fill="gains")) +
  geom_area(aes(y=-intermediate_loss, alpha="intermediate", fill="losses")) + 
  geom_hline(yintercept=0) + 
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers), 
             aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(.~chr, scales="free_x", space="free_x") + 
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) + 
  scale_alpha_manual(values=c("clonal"=1, "intermediate"=0.5, "tip"=0.2)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  scale_x_continuous(breaks=c(0), labels="") + 
  scale_y_continuous(breaks=seq(-ylm, ylm, by=20), 
                     labels=abs(seq(-ylm, ylm, by=20)), 
                     limits=c(-ylm,ylm)) + 
  xlab("") + ylab("Number of cases") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}

# Save the output
ggsave("results/4_clonality_frequency/genome_wide_clonality_intermediate_tree_based_frequency.png", plot_cna_recurence, width=11, height=2.1)

ylm = 60

plot_cna_recurence=
  data_per_window_summary %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)]) %>%
  mutate(chr=factor(chr, c(1:22,"X","Y"), ordered=TRUE)) %>%
  ggplot(aes(x=pos, y=tip_gain)) + 
  geom_area(aes(alpha="tip", fill="gains")) +
  geom_area(aes(y=-tip_loss, alpha="tip", fill="losses")) + 
  geom_hline(yintercept=0) + 
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers), 
             aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(.~chr, scales="free_x", space="free_x") + 
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) + 
  scale_alpha_manual(values=c("clonal"=1, "intermediate"=0.5, "tip"=0.2)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  scale_x_continuous(breaks=c(0), labels="") + 
  scale_y_continuous(breaks=seq(-ylm, ylm, by=20), 
                     labels=abs(seq(-ylm, ylm, by=20)), 
                     limits=c(-ylm,ylm)) + 
  xlab("Chromosomes") + ylab("Number of cases") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}

# Save the output
ggsave("results/4_clonality_frequency/genome_wide_clonality_tip_tree_based_frequency.png", plot_cna_recurence, width=11, height=2.1)
