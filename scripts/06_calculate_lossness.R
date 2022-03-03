# Script to calculate lossness metrics

# Functions
source("scripts/00_general_functions.R")

# Where is the low pass data?
pipeline_root = "~/remotes/forecast_low_pass/"

# Minimum PGA
min_pga = 0.01

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Get all files
call_files = list.files("results/0_min_het_best_fit_ploidy_search/calls/", 
                        pattern = "_cna_ploidy_search_calls.txt",
                        full.names = T)

# Find the call files
per_patient_samples = lapply(patients, function(p) call_files[grep(p, x = call_files)])

# Read in the data
per_patient_calls = lapply(per_patient_samples, FORECAST_read.files)

# Arm boundaries
armb = read.table("refs/chrArmBoundaries_hg38.txt", header = T, stringsAsFactors = F)

# Bin names
bin_names_all = extractChrRows(rownames(per_patient_calls[[1]][[1]]))

# Run through chromosomes and get telo and centro bins
bad_pos = unlist(lapply(c(1:22,"X","Y"), function(chr) {
  
  start_bin = paste0(chr,"_",min(bin_names_all[bin_names_all$chr==chr,"start"]))
  end_bin   = paste0(chr,"_",max(bin_names_all[bin_names_all$chr==chr,"end"]))
  
  mid_loc   = armb[armb$Chr==chr,2]
  
  q_locs = bin_names_all[bin_names_all$chr==chr,"start"][bin_names_all[bin_names_all$chr==chr,"start"] - mid_loc > 0]
  p_locs = bin_names_all[bin_names_all$chr==chr,"end"][bin_names_all[bin_names_all$chr==chr,"end"] - mid_loc < 0]
  
  start_q_arm = paste0(chr,"_",min(q_locs))
  end_p_arm   = ifelse(length(p_locs)>0, paste0(chr,"_",max(p_locs)), "")
  
  res = c(start_bin, end_p_arm, start_q_arm, end_bin)
  
}))

# Read in the patient ploidies
patient_ploidies = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                              header = T, stringsAsFactors = F)
patient_ploidies$Ploidy = round(patient_ploidies$PsiT)

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", 
                             header = T, stringsAsFactors = F)

# Make a blacklist of impure samples
normal_samples = per_sample_mets$PGA < min_pga | per_sample_mets$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", per_sample_mets$Sample) | metrics$Sample=="IM1017_B2_DW1"
normal_samples = per_sample_mets$Sample[normal_samples]

# Get indices of the tumour samples
per_patient_calls = lapply(per_patient_calls, function(p) {
  
  t_ind = which(unlist(lapply(p, function(n) {!removeBamNameEnding(colnames(n)) %in% normal_samples})))
  
  return(p[t_ind])
  
})

# Now some patients may have been removed due to having no CNAs
patient_no_data   = which(unlist(lapply(per_patient_calls, length))>0)
per_patient_calls = per_patient_calls[patient_no_data]
patients          = patients[patient_no_data]

# Get the chr per bin vector
chrs = unlist(lapply(strsplit(rownames(per_patient_calls[[1]][[1]]), split = ":"), function(i) i[1]))

# % calls that are losses
pplossness = lapply(1:length(per_patient_calls), function(i) {
  
  pcs = per_patient_calls[[i]]
  pat = patients[i]
  
  bp = patient_ploidies[patient_ploidies$Case==pat,"Ploidy"]
  
  bp = ifelse(chrs%in%c("X","Y"), 
              bp/2, 
              bp)
  
  print(tail(bp))
  
  res = lapply(pcs, function(s) fraction(s[,1]<bp))
  
  names(res) = unlist(lapply(per_patient_calls[[i]], colnames))
  
  return(res)
  
})

# Mean loss %
mlossness  = lapply(pplossness, function(p) mean(unlist(p)))
names(mlossness) = patients

# Collapse segments
per_patient_collapsed_segs = lapply(per_patient_calls, function(p) lapply(p, collapseAsSegments, chrs = c(1:22,"X","Y")))
names(per_patient_collapsed_segs) = patients

# % segments that are losses
per_patient_cs_lossness = lapply(1:length(per_patient_collapsed_segs), function(i) {
  
  p = per_patient_collapsed_segs[[i]]
  pat = patients[i]
  
  bp = patient_ploidies[patient_ploidies$Case==pat,"Ploidy"]
  
  res = unlist(lapply(p, function(s) {
    
    bp = ifelse(s$chr%in%c("X","Y"), 
                bp/2, 
                bp)
        
    res = fraction(s$CN<bp & !paste0(s$chr,"_",s$start) %in% bad_pos & !paste0(s$chr,"_",s$end) %in% bad_pos)
    
    return(res)
    
  }))
  
  names(res) = unlist(lapply(per_patient_calls[[i]], colnames))
  
  return(res)
  
})

# Add names
names(per_patient_cs_lossness) = patients

# Write it out as an RDS file
saveRDS(per_patient_cs_lossness, file = "results/6_lossness/FORECAST_collapsed_segment_lossness.rds")

# Mean loss %
mcslossness  = lapply(per_patient_cs_lossness, function(p) mean(unlist(p)))
names(mcslossness) = patients

# Make a plot too
pdf("results/6_lossness/FORECAST_per_patient_lossness.pdf", width = 21)
boxplot(per_patient_cs_lossness[order(unlist(mcslossness), decreasing = T)], 
        las = 2, pch = 20, cex = 0.5, ylab = "% segments loss", main = "FORECAST Lossness Metric - Searching for BRCAness")
dev.off()
