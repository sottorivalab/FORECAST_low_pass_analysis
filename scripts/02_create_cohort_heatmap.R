# A script to read in the calls and create a matrix and then a heatmap
library(ComplexHeatmap)
library(circlize)

# Functions
source("scripts/00_general_functions.R")

# Where is the low pass data?
pipeline_root = "~/remotes/forecast_low_pass/"

# Minimum PGA to use
min_pga = 0.01
# Include X,Y?
inc_sex = T

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
normal_samples = per_sample_mets$PGA < min_pga | per_sample_mets$Assessment=="Failed" | grepl("_N[0-9]{1,2}_", per_sample_mets$Sample) | per_sample_mets$Sample=="IM1017_B2_DW1"
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

# Convert to a matrix
per_patient_call_mat = lapply(per_patient_calls, function(p) {
  
  m = do.call(cbind, p)
  
  colnames(m) = gsub("_DW1_dups_valid", "", colnames(m))
  
  return(m)
  
})

# Collapse lists into a matrix
per_patient_call_mat = lapply(per_patient_call_mat, function(p) {
  
  if(ncol(p)>1) {
  
    x = p[,hclust(dist(t(p)))$order]
    
  } else {x = p}
  
  return(x)
  
})

# Make a matrix of the whole cohort
mat = do.call(cbind, per_patient_call_mat)

# Make a patient vector
pats = extractPatientName(colnames(mat))

#Patient separations
pat.ends   = c(0, cumsum(unlist(lapply(unique(pats), function(x) count(pats == x)))))

# Extract chromosome vector
chrs = unlist(lapply(strsplit(rownames(mat), split = "[:]"), function(i) i[1]))

#Chromosome separation positions
clp_chr  = ifelse(inc_sex, 24, 22)
chr.ends = cumsum(rle(chrs)$lengths)[-clp_chr]
if(inc_sex) {lab_chr = c(1:22,"X","Y")} else {lab_chr = c(1:22)}

# Annotate bottom with chromosomes
bot.column.anno = columnAnnotation(link = anno_mark(at = round(c(0, chr.ends) + rle(chrs)$lengths / 2), 
                                                    labels = lab_chr, 
                                                    side = "bottom",
                                                    link_width = unit(0.25, "cm")))
# Find the middle
pannot = round((table(pats) / 2) + c(0, cumsum(table(pats))[-length(table(pats))]))

# Add patient naming
left.row_anno = rowAnnotation(link = anno_mark(at = pannot,
                                               labels = names(pannot),
                                               side = "left",
                                               link_width = unit(0.5, "cm")),
                              width = unit(3, "cm"))

# Output the heatmap
pdf("results/1_heatmap/Low_pass_cohort_heatmap.pdf", height = 20, width = 12)

# plot the heatmap
Heatmap(t(mat), cluster_rows = F, cluster_columns = F, show_column_names = F, 
        bottom_annotation = bot.column.anno,
        left_annotation = left.row_anno,
        show_row_names = F, col = colorRamp2(c(0, 2, 4, 8, 16), 
                                             c("blue", "white", "red", "darkred", "black")), name = "CN")

#Add lines
for(boundary in pat.ends / length(pats)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(0, 1), c(1 - boundary, 1 - boundary), gp = gpar(lwd = 0.5))
  })
  
}

#Add lines
for(boundary in c(0,1)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lwd = 1))
  })
  
}

#Add lines
for(boundary in chr.ends / nrow(mat)) {
  
  #Add the lines
  decorate_heatmap_body("CN", {
    grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
  })
  
}

dev.off()
