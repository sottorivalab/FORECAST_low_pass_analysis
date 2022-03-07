# Code to make supplementary trees with gleason grade annotated and mutations listed

library(ggtree)
library(ggtext)
library(ggpubr)
library(circlize)

# Targeted data repo
targeted_repo = "target_repo/"

# Roots for data input
medicc_root     = "~/remotes/misc/"

# Colour scale for gleason
gleason_cols = colorRamp2(c(1,2,3,4,5), colors = c("#009E73", "#F0E442", "#D55E00", "#5D3A9B", "#000000"))

# Read in the table of mutations contained in this patient
mutations = read.table(paste0(targeted_repo,"results/mutation_calling/Per_sample_gene_mutation_status_with_change.txt"), 
                       sep = "\t", header = T, stringsAsFactors = F)

# Add in the indels too!
indels = read.table(paste0(targeted_repo,"results/mutation_calling/Per_sample_gene_mutation_status_with_change_indels.txt"),
                    sep = "\t", header = T, stringsAsFactors = F)

# Combine them
mutations = rbind(mutations, indels)

# Get the samples sequenced
sequenced_samples = read.table(paste0(targeted_repo,"results/mutation_calling/Per_sample_gene_mutation_status.txt"),
                               sep = "\t", header = T, stringsAsFactors = F)
sequenced_samples = sequenced_samples[!is.na(sequenced_samples$AKT1),]
seq_samples       = paste0(sequenced_samples$Patient,"_",sequenced_samples$Sample)

# Get per region gleason scores
gleason_cont = read.csv("results/8_collect_metrics/Gleason_per_sample_digital_pathology.csv")

# Get patients
medicc_patients = gsub("_regions.tsv", "", list.files(paste0(medicc_root,"forecast"), pattern = "_regions.tsv"))

# Run through them
for(patient in medicc_patients) {
  
  # Get patient gleason
  p_gleason = gleason_cont[gleason_cont$Patient==patient,]
  
  # Read in the file
  tree = read.tree(paste0(medicc_root,"forecast/",patient,"_regions_final_tree.new"))
  
  # Make a short df of tree samples
  sample_gleason = data.frame(Sample = tree$tip.label)
  
  # Get cont_gleason
  sample_gleason = merge(sample_gleason, p_gleason[,c("Sample", "Grade.Group", "Cont_Gleason")], all.x = T, sort = F)
  
  # Get the colour
  sample_gleason$colour = gleason_cols(sample_gleason$Grade.Group)
  
  # Make NAs some colour
  sample_gleason$colour[is.na(sample_gleason$colour)] = "grey"
  sample_gleason$colour[sample_gleason$Sample=="diploid"] = "#009E73"
  
  # Replace the patient and DW1
  tree$tip.label = gsub("_DW1", "", gsub(paste0(patient,"_"), "", tree$tip.label))
  
  # Basic tree structure
  p = ggtree(tree, size = 1) + geom_tiplab(hjust = -0.5, size = 8) + 
    geom_tippoint(col = sample_gleason$colour, size = 6) +
    ggtitle(paste0(patient," (CNA)")) + theme(title=element_text(size=20, colour = "gray20"), plot.title = element_text(hjust = 0.5))
  
  # Add events labels?
  p = p + geom_label(data = p$data[p$data$branch.length > 0,], aes(x = branch, label = branch.length), label.size = 0, nudge_y = 0.2)
  
  # Customise width?
  p2 = p
  p  = p2 + ggplot2::xlim(0, max(ape::node.depth.edgelength(tree)) * 1.1)
  
  # Get patient mutations
  p_muts = mutations[mutations$Patient==patient,]
  
  # Sequenced samples
  p_seq_samples = seq_samples[grep(paste0(patient,"_"), seq_samples)]
  p_seq_samples = gsub(paste0(patient,"_"), "", p_seq_samples)
  
  # Remove those with VAF < 0.05 
  p_muts = p_muts[p_muts$VAF >= 0.05,]
  
  # Run through mutations and gather which samples have the mutation
  muts = unlist(lapply(sort(unique(p_muts$Protein)), function(i) {
    
    res = paste0(gsub("_DW1", "", sort(p_muts[p_muts$Protein==i,"Sample"])), collapse = ", ")
    
    paste0(c(i, res), collapse = ": ")
    
  }))
  
  # Get things in the right format
  cna_not_seq = p$data$label[p$data$isTip][!p$data$label[p$data$isTip] %in% gsub("_DW1", "", p_seq_samples)]
  
  cna_not_seq = cna_not_seq[cna_not_seq!="diploid"]
  
  cna_not_seq = paste0(cna_not_seq, collapse = ", ")
  
  not_in_cna_tree = gsub("_DW1", "", p_seq_samples)[!gsub("_DW1", "", p_seq_samples) %in% p$data$label[p$data$isTip]]
  
  not_in_cna_tree = not_in_cna_tree[not_in_cna_tree!="BC1_DNA1"]
  
  not_in_cna_tree = paste0(not_in_cna_tree, collapse = ", ")
  
  if(cna_not_seq!="") {muts = c(muts, paste0(cna_not_seq," TES unavailable"))}
  if(not_in_cna_tree!="") {muts = c(muts, paste0(not_in_cna_tree," not included in tree"))}
  
  muts = paste0(muts, collapse = "<br><br>")
  
  muts = gsub("[*]", "\\\\*", muts)
  
  df = data.frame(
    label = muts,
    x = c(0),
    y = c(0.8),
    hjust = c(0),
    vjust = c(1),
    orientation = c("upright"),
    color = c("black"),
    fill = c("white")
  )
  
  p2 = ggplot(df) +
    aes(
      x, y, label = label, color = color, fill = fill,
      hjust = hjust, vjust = vjust,
      orientation = orientation
    ) +
    geom_textbox(width = unit(0.8, "npc"), size = 6) +
    scale_discrete_identity(aesthetics = c("color", "fill", "orientation")) +
    xlim(0, 1) + ylim(0, 1) + xlab("") + ylab("") + theme_void()
  
  p_res = ggarrange(p, p2,
            ncol = 2, nrow = 1)
  ggsave(plot = p_res, paste0("results/7_medicc_analysis/",patient,"_medicc_tree_mutation_legend.png"), width = 16, height = 6)

}

