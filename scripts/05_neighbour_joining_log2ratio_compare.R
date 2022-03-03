# Assess purity adjusted log2ratio similarity and then make NJ trees
library(phangorn)
library(phytools)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

# Functions
source("scripts/00_general_functions.R")

# Where is the low pass data?
pipeline_root = "~/remotes/forecast_low_pass/"

# What is the minimum required PGA to be tumour?
min_pga = 0.01

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Find the seg files
per_patient_samples   = lapply(patients, FORECAST_list.files, pattern = "_cna_segments.txt")

# Read in the data
per_patient_segs = lapply(per_patient_samples, FORECAST_read.files)

# Find the bin files
per_patient_samples   = lapply(patients, FORECAST_list.files, pattern = "_bins.txt")

# Read in the data
per_patient_bins = lapply(per_patient_samples, FORECAST_read.files)

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", header = T, stringsAsFactors = F)

# Make a blacklist of impure samples
normal_samples = per_sample_mets$PGA < min_pga | per_sample_mets$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", per_sample_mets$Sample) | per_sample_mets$Sample=="IM1017_B2_DW1"
normal_samples = per_sample_mets$Sample[normal_samples]

# Plot the distribution of samples per patient
hist(unlist(lapply(per_patient_samples, length)), breaks = 100, xlim = c(0,25), 
     xlab = "Samples", main = "Distribution of samples per patient")

# Get indices of the tumour samples
per_patient_segs = lapply(per_patient_segs, function(p) {
  
  t_ind = which(unlist(lapply(p, function(n) {!removeBamNameEnding(colnames(n)) %in% normal_samples})))
  
  return(p[t_ind])
  
})

# Get indices of the tumour samples
per_patient_bins = lapply(per_patient_bins, function(p) {
  
  t_ind = which(unlist(lapply(p, function(n) {!removeBamNameEnding(colnames(n)[5]) %in% normal_samples})))
  
  return(p[t_ind])
  
})

# Number of samples per patient after removal of low PGA samples
hist(unlist(lapply(per_patient_segs, length)), breaks = 100, xlim = c(0,25), 
     xlab = "Samples", main = "Distribution of samples per patient (tumour)")

# Now some patients may have been removed due to having no CNAs
patient_no_data  = which(unlist(lapply(per_patient_segs, length))>0)
per_patient_segs = per_patient_segs[patient_no_data]
per_patient_bins = per_patient_bins[patient_no_data]
patients         = patients[patient_no_data]

# Calculate spearman coefficients
per_patient_cor_mat = lapply(per_patient_bins, function(p) {
  
  if(length(p) > 1) {
  
    combinations = combn(length(p), 2)
    
    comb_dists = apply(combinations, 2, function(i) {
      
      cor(p[[i[1]]][,5], p[[i[2]]][,5], method = "spearman")
      
    })
    
    mat = matrix(nrow = length(p), ncol = length(p))
    
    for(i in 1:ncol(combinations)) {
      
      mat[combinations[2,i],combinations[1,i]] = comb_dists[i]
      
    }
    
    samples = unlist(lapply(p, function(s) removeBamNameEnding(colnames(s)[5])))
    
    patient = extractPatientName(samples[1])
    
    samples = gsub(paste0(patient,"_"), "", samples)
    
    samples = gsub("_DW1", "", samples)
    
    rownames(mat) = samples
    colnames(mat) = samples
    
    return(mat)
    
  } else {return(NULL)}
  
})

# Melt it
per_patient_pairwise = lapply(1:length(patients), function(i) {
  
  res = per_patient_cor_mat[[i]]
  
  if(!is.null(res)) {
    
    colnames(res) = paste0(patients[i],"_",colnames(res),"_DW1")
    rownames(res) = paste0(patients[i],"_",rownames(res),"_DW1")
    
    res = melt(res)
    
    res = na.omit(res)
    
    res$Patient = patients[i]
    
  } else {res = NULL}
  
  return(res)
  
})

# Combine
per_patient_pairwise = do.call(rbind, per_patient_pairwise)

# Name columns
colnames(per_patient_pairwise)[1:3] = c("S1", "S2", "Spearman")

# Output the pair data
write.csv(per_patient_pairwise, file = "results/5_boxplot_log2ratio_dist/Spearman_rho_per_patient_pairwise.csv", 
          quote = F, row.names = F)

# Collapse
per_patient_cors        = lapply(per_patient_cor_mat, function(p) na.omit(as.vector(p)))
names(per_patient_cors) = patients

# Write out the Spearman distance information
saveRDS(per_patient_cors, file = "results/5_boxplot_log2ratio_dist/Spearman_rho_per_patient.rds")

# Add on a false normal
per_patient_segs = lapply(per_patient_segs, function(p) {
  
  p[[length(p)+1]] = cbind(rep(0, times = nrow(p[[1]])))
  
  colnames(p[[length(p)]]) = "Normal"
  
  return(p)
  
})

# Apply distance measure to all pairs in sample
per_patient_dist_mat = lapply(per_patient_segs, function(p) {
  
  combinations = combn(length(p), 2)
  
  comb_dists = apply(combinations, 2, function(i) {
    
    log2ratio_comparison(segs_col_a = p[[i[1]]][,1], 
                         segs_col_b = p[[i[2]]][,1], 
                         normalise_to_exp = T, exp_distance = 200)
    
  })
  
  mat = matrix(nrow = length(p), ncol = length(p))
  
  for(i in 1:ncol(combinations)) {
    
    mat[combinations[2,i],combinations[1,i]] = comb_dists[i]
    
  }
  
  samples = unlist(lapply(p, function(s) removeBamNameEnding(colnames(s))))
  
  patient = extractPatientName(samples[1])
  
  samples = gsub(paste0(patient,"_"), "", samples)
  
  samples = gsub("_DW1", "", samples)
  
  rownames(mat) = samples
  colnames(mat) = samples
  
  return(mat)
  
})

# Removing the normal comparison
ds = lapply(per_patient_dist_mat, function(p) p[-nrow(p),-ncol(p)])

# Extract distances as a vector
ds = lapply(ds, function(d) na.omit(as.vector(d)))

# What's the order
dorder = order(unlist(lapply(ds, mean)), decreasing = T)
# Make a boxplot of the intratumour distances
pdf("results/5_boxplot_log2ratio_dist/FORECAST_cohort_log2ratio_similarities.pdf", width = 21)
boxplot(ds[dorder], names = patients[dorder], las = 2, pch = ".", ylab = "L2RSS")
dev.off()

# Add names
names(ds) = patients

# Write out the log2ratio distance information
saveRDS(ds, file = "results/5_boxplot_log2ratio_dist/L2RSS_scores_per_patient.rds")

# Make a boxplot of the intratumour distances
pdf("results/5_boxplot_log2ratio_dist/FORECAST_cohort_spearman_similarities.pdf", width = 21)
boxplot(per_patient_cors[dorder], names = patients[dorder], las = 2, pch = ".", ylab = "Spearman", ylim = c(0,1))
dev.off()
