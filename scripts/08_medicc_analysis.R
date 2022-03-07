# This script analyses the MEDICC2 output

library(phytools)
library(usedist)
library(ggplot2)
library(cowplot)
library(adephylo)
library(ComplexHeatmap)
library(circlize)
options(scipen = 999)

# Roots for data input and output
medicc_root     = "~/remotes/misc/"
pipeline_root   = "~/remotes/forecast_low_pass/"

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Loop them
trees = lapply(patients, function(patient) {
  
  tree = try(read.tree(paste0(medicc_root,"forecast/",patient,"_regions_final_tree.new")))
  
  return(tree)
  
})

# Loop them - I edited MEDICC2 to output this
mrcas = lapply(patients, function(patient) {
  
  mrca = try(read.table(paste0(medicc_root,"forecast/",patient,"_regions_mrca.txt"), 
                        stringsAsFactors = F, header = F)$V1)
  
  if(class(mrca)=="try-error") {mrca = NA}
  
  return(mrca)
  
})

# Loop them
total_events = lapply(trees, function(tree) {
  
  if(class(tree)=="phylo") {events = sum(tree$edge.length)} else {events = NA}
  
  return(events)
  
})

# Loop them
events_breakdown = lapply(trees, function(tree) {
  
  if(class(tree)=="phylo") {breakdown = data.frame(A = tree$edge[,1],
                                                   B = tree$edge[,2],
                                                   Edge_Length = tree$edge.length, 
                                                   A_Label = c(tree$tip.label, tree$node.label)[tree$edge[,1]], 
                                                   B_Label = c(tree$tip.label, tree$node.label)[tree$edge[,2]])
  
  } else {breakdown = NULL}
  
  return(breakdown)
  
})

# Collect tip lengths
tip_lengths = do.call(rbind, lapply(events_breakdown, function(i) i[,c("Edge_Length", "B_Label")]))

# Filter for things we don't need
tip_lengths = tip_lengths[!grepl("internal_", tip_lengths$B_Label),]
tip_lengths = tip_lengths[!grepl("diploid", tip_lengths$B_Label),]

# Format
tip_lengths = tip_lengths[,2:1]
colnames(tip_lengths) = c("Sample", "Tip_Length")

# Output that
write.table(tip_lengths, file = "results/7_medicc_analysis/Tip_lengths_per_sample_medicc2.csv", quote = F, sep = ",", row.names = F)

# Get clonal events
clonal_events = lapply(events_breakdown, function(i) {
  
  ce = i[i$A_Label=="" & i$B_Label!="diploid","Edge_Length"]
  
  if(length(ce)==0) {ce = NA}
  
  return(ce)
  
})

# Loop them
event_pairs = lapply(trees, function(tree) {

  if(class(tree)=="phylo") {cna_dists = distTips(tree)} else {cna_dists = NA}
  
  return(cna_dists)
  
})

# Loop them
summary = lapply(patients, function(patient) {
  
  summary = try(read.table(paste0(medicc_root,"forecast/",patient,"_regions_summary.tsv"), 
                           stringsAsFactors = F, header = F, sep = "\t", fill = T))
  
  if(class(summary)=="try-error") {summary = NA}
  
  return(summary)
  
})

# Loop them
region_status = lapply(patients, function(patient) {
  
  cnas = try(read.table(paste0(medicc_root,"forecast/",patient,"_regions_final_cn_profiles.tsv"), 
                        stringsAsFactors = F, header = T, sep = "\t", fill = T))
  
  if(class(cnas)=="try-error") {cnas = NA}
  
  return(cnas)
  
})

# Read in the per tree summaries and export it
tree_summaries = data.frame(Patient = patients,
                            MRCA = unlist(mrcas),
                            Total_Events = unlist(total_events),
                            Clonal_Events = unlist(clonal_events),
                            Subclonality = (unlist(total_events) - unlist(clonal_events)) / unlist(total_events))

# Make a ggplot
ggplot(data = tree_summaries, aes(x = Subclonality)) + 
  geom_histogram(bins = round(diff(range(tree_summaries$Subclonality, na.rm = T)) * 100), color = "black", fill = "gray") + 
  geom_vline(aes(xintercept = mean(Subclonality, na.rm = T)), linetype = "dashed", size = 0.6) + 
  geom_vline(aes(xintercept = quantile(Subclonality, probs = 1/3, na.rm = T)), linetype = "dotted", size = 0.6) +
  ggtitle(paste0("Subclonality of all samples in the cohort (mean = ",round(mean(tree_summaries$Subclonality, na.rm = T), digits = 2),")")) +
  ylab("Count") +
  theme_cowplot() + theme(plot.title = element_text(hjust = 0.5))
ggsave("results/7_medicc_analysis/Histogram_all_subclonality_cohort.pdf", height = 6, width = 7)

# Write it out
write.table(tree_summaries, file = "results/7_medicc_analysis/medicc_analysis_summary.txt", 
            sep = "\t",
            row.names = F, quote = F)

# Read in the pipeline metrics
metrics    = read.table("results/2_metrics/metric_assessment_table.txt", header = T, stringsAsFactors = F)

# Merge with metrics the tip lengths to look for a correlation there
metrics_tip = merge(metrics, tip_lengths, by = "Sample")

# Plot it to show outliers
ggplot(metrics_tip, aes(x = Purity, y = Tip_Length)) + geom_point() + theme_cowplot() + 
  geom_smooth() + theme(plot.title = element_text(hjust = 0.5))
