# Script to make input for MEDICC2

library(reshape2)

# Functions
source("scripts/00_general_functions.R")

# Roots for data input and output
medicc_root     = "~/remotes/misc/"
pipeline_root   = "~/remotes/forecast_low_pass/"

# Minimum proportion genome altered allowed
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

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", header = T, stringsAsFactors = F)

# Make a blacklist of samples
normal_samples = per_sample_mets$PGA < min_pga | per_sample_mets$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", per_sample_mets$Sample) | metrics$Sample=="IM1017_B2_DW1"
normal_samples = per_sample_mets$Sample[normal_samples]

# Get indices of the tumour samples
per_patient_calls = lapply(per_patient_calls, function(p) {
  
  t_ind = which(unlist(lapply(p, function(n) {!removeBamNameEnding(colnames(n)) %in% normal_samples})))
  
  return(p[t_ind])
  
})

# Make it a matrix
per_patient_mats        = lapply(per_patient_calls, function(p) do.call(cbind, p))
names(per_patient_mats) = patients

# Count the number of samples
n_samples = unlist(lapply(per_patient_mats, ncol))

# Run though the patients
for(p in patients) {
  
  # Can only do it with multiple samples
  if(n_samples[p]>1) {
  
    # Make MEDICC2 output
    m = per_patient_mats[[p]]
    
    # Split out row names by chr, start, end
    locations_split = lapply(strsplit(rownames(m), "-"), function(i) unlist(strsplit(i, split = ":")))
    
    # Add location columns
    chr   = unlist(lapply(locations_split, function(i) i[1]))
    start = as.numeric(unlist(lapply(locations_split, function(i) i[2])))
    end   = as.numeric(unlist(lapply(locations_split, function(i) i[3])))
    
    # Make a dataframe
    chr_pos = data.frame(chr = chr,
                         start = start,
                         end = end)
    
    # Make an identifier
    ids = sapply(1:nrow(m), function(i) paste0(c(chr[i], m[i,]), collapse = "_"))
    
    # Collapse
    run_len = rle(ids)
    
    # Cum summed
    cumulative_len = cumsum(run_len$lengths)
    
    # Create a matrix
    res = data.frame(chrom = chr[cumulative_len],
                     start = start[c(1, cumulative_len[-length(cumulative_len)]+1)],
                     end = end[cumulative_len],
                     m[cumulative_len,])
    
    # For diploid root
    res_root = res[,1:3]
    res_root$cn = ifelse(res_root$chrom %in% c("X", "Y"), 1, 2)
    res_root$sample_id = "diploid"
    res_root = res_root[,c(5,1:3,4)]
    rownames(res_root) = NULL
    
    # Rearrange
    res = cbind(res[,1:3], melt(res[,4:ncol(res)]))
    res = res[,c(4,1:3,5)]
    colnames(res) = c("sample_id", "chrom", "start", "end", "cn")
    
    options(scipen = 999)
    
    res = rbind(res, res_root)
    
    # Write it out
    write.table(res, file = paste0(medicc_root,"forecast/",p,"_regions.tsv"),
                sep = "\t", row.names = F, quote = F)
    
  }
  
}
