# Functions

# Range function for consensus ploidy analysis
range01 = function(x){(x-min(x))/(max(x)-min(x))}

# Remove zero function
removeZero = function(x) {x[x!=0]}

# Get tertiles
tertile = function(x) {quantile(x, probs = c(1/3, 2/3))}

# Remove a column's NAs
colNArm = function(x, sub) {x[!is.na(x[,sub]),]}

# Convert DP naming
convertPathSampleName = function(x) {
  
  # Edit patient/sample names
  x = gsub(".csv", "", x)
  x = unlist(lapply(strsplit(x, split = ""), function(i) paste0(c(i[1:6],"_",i[7:length(i)],"_DW1"), collapse = "")))
  
  return(x)
  
}

# Make a summary list
makeSummaryList = function(x) {
  
  res = list()
  length(res) = length(x)
  names(res) = x
  
  return(res)
  
}

# Check which patients have a mutation applying NAs where needed
applyMutStatus = function(query_list, mutated_cases, attempted_seq, bad_cases) {
  
  res = query_list %in% mutated_cases
  res[!query_list %in% attempted_seq] = NA
  res[query_list %in% bad_cases] = NA
  
  return(res)
  
}

# Make location dataframe from chr pos rownames
extractChrRows = function(x) {
  
  m = do.call(rbind, strsplit(x, split = "[:]|-"))
  
  df = data.frame(chr = m[,1],
                  start = as.integer(m[,2]),
                  end = as.integer(m[,3]))
  
  df$chr = factor(df$chr, levels = c(1:22, "X", "Y"))
  
  return(df)
  
}

# Remove generic bam name end
removeBamNameEnding = function(x) {gsub("_dups_valid", "", x)}

# Extract patient name
extractPatientName = function(x) {unlist(lapply(strsplit(x, split = "_"), function(i) i[1]))}

# A sum of squared difference comparison of segmented log2ratio per bin, normalising for purity differences 
log2ratio_comparison = function(segs_col_a, segs_col_b, exp_distance = 1848.691, normalise_to_exp = T, min_purity = 0.2) {
  
  # Calculate contineous copy number
  calcCN = function(lrrs, rho, psit, gamma = 1) {
    
    psi = (2*(1 - rho)) + (rho*psit)
    
    n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
    
    return(n)
    
  }
  
  # What is our parameter search of purities?
  parameter_comparison = rbind(cbind(seq(min_purity, 0.99, by = 0.01), 1),
                               cbind(1, seq(1, min_purity, by = -0.01)))
  
  # Here we do a search of purity pairs
  search = lapply(1:nrow(parameter_comparison), function(r) {
    
    # Selected parameters for iteration
    rhoA = parameter_comparison[r,1]
    rhoB = parameter_comparison[r,2]
    
    # Continuous copy number calculation 
    CNa  = calcCN(lrrs = segs_col_a, rho = rhoA, psit = 2)
    CNb  = calcCN(lrrs = segs_col_b, rho = rhoB, psit = 2)
    
    # Sum of squared differences (maybe normalise for number of bins?)
    dist = sum((CNa - CNb)^2)
    
    return(dist)
    
  })
  
  # Distance results for parameter comparisons
  res = cbind(parameter_comparison, unlist(search))
  
  # Which has the shortest distance
  R = which.min(res[,3])
  
  # Get the d
  d = res[R,3]
  
  if(normalise_to_exp) {
    
    # Normalise the distance to the cohort (hard coded for now)
    d = d / exp_distance # This number is the median dist in non-same patient comparisons
    if(d>1) {d = 1} # Cap at 1
    
  }
  
  return(d)
  
}

# List files in the FORECAST structure
FORECAST_list.files = function(p, pattern, sub_dir = "/QDNASEQ") {
  
  files = list.files(paste0(davros_mount,"data/",p,sub_dir), 
                     recursive = T, pattern = pattern, full.names = T)
  
  return(files)
  
}

# Make a function for reading in files
FORECAST_read.files = function(p) {
  
  dat = lapply(p, function(f) {
    
    try(read.table(f, header = T, stringsAsFactors = F))
    
  })
  
}

# Speed up finding files in the lpWGS pipeline
findFileslpWGS = function(file_pattern, sub_dir, pipeline_root = "~/remotes/forecast_low_pass/data/", patients = "auto") {
  
  if(patients == "auto") {
    
    dirs     = list.files(pipeline_root)
    
    patients = c(dirs[grep("^DNT", dirs)], dirs[grep("^IM", dirs)])
    
    print(patients)
    
  }
  
  per_patient = lapply(patients, function(p) {
    
    list.files(paste0(pipeline_root,p,"/",sub_dir), 
               pattern = file_pattern, recursive = T, full.names = T)
    
  })
  
  files = unlist(per_patient)
  
  return(files)
  
}

# Save here a function for calculating inter-patient divergences to use as baseline divergence
FORECAST_L2RSS_exp_calc = function(its = 1000) {
  
  unlist(lapply(1:its, function(i) {
    
    cases = sample(1:length(per_patient_segs), size = 2)
    
    p1 = per_patient_segs[[cases[1]]]
    p2 = per_patient_segs[[cases[2]]]
    
    p1 = p1[[sample(1:length(p1), size = 1)]]
    p2 = p2[[sample(1:length(p2), size = 1)]]
    
    d = log2ratio_comparison(p1[,1], p2[,1], normalise_to_exp = F)
    
    return(d)
    
  }))
  
}

# Get p and q arm assignment
getPQ = function(df, arm_pos) {
  
  psnqs = unlist(lapply(1:nrow(df), function(i) {
    
    r = df[i,]
    
    pos = (r$end - r$start) + r$start
    
    chr = r$chromosome
    
    b = arm_pos[which(arm_pos[,1] == chr),2]
    
    if(pos <  b) {a = "p"}
    if(pos >= b) {a = "q"}
    
    return(a)
    
  }))
  
  return(psnqs)
  
}

# Summarise list to dataframe
list_summary_df = function(l, value_name, summary = c("mean", "median", "sum")) {
  
  summary = match.arg(summary)
  
  if(summary=="mean")   {collapsed_data = unlist(lapply(l, mean))}
  if(summary=="median") {collapsed_data = unlist(lapply(l, median))}
  if(summary=="sum")    {collapsed_data = unlist(lapply(l, sum))}
  
  df = data.frame(Patient = names(collapsed_data), value = collapsed_data)
  
  colnames(df)[2] = value_name
  
  rownames(df) = NULL
  
  return(df)
  
}

# Collapse calls to segments
collapseAsSegments = function(calls, chrs = 1:22) {
  
  do.call(rbind, lapply(chrs, function(chr) {
    
    index = grep(paste0("^", chr,":"), rownames(calls), perl = T)
    
    bin_names = rownames(calls)[index]
    
    chr_cn = calls[index,1]
    
    rle_res = rle(chr_cn)
    
    loc_annot = extractChrRows(bin_names)
    
    data.frame(chr = chr, start = loc_annot$start[c(1, (cumsum(rle_res$lengths)+1)[-length(rle_res$lengths)])], 
               end = loc_annot$end[cumsum(rle_res$lengths)], CN = rle_res$values)
    
  }))
  
}

# Make a function called count
count = function(x) {length(which(x))}

# Make a function that just does a fraciton
fraction = function(x) {length(which(x)) / length(x)}

# This function normalise expected signal from male sex chromosomes
male_normalise_lrr = function(logr, chr_vec) {log(ifelse(chr_vec %in% c("X","Y"), 2, 1) * 2^logr, base = 2)}

# Call on X chromosome given purity and ploidy known already
callXchromsome = function(sex_lrrs, psi, psit, purity) {

  # Adjust them so that they are normalised to what is expected
  hap_adj_lrrs = log(((2^sex_lrrs) * psi) / (psi/2), base = 2)
  
  # This is the haploid expected state, i.e. 1
  hap_e = 1
  
  # Half the PsiT to represent haploidy
  hap_psit = psit / 2
  
  # Get psi
  hap_psi = (hap_e*(1 - purity)) + (purity*hap_psit)
  
  # Calculate using purity
  n = ((hap_psi*(2^(hap_adj_lrrs/1))) - (hap_e * (1 - purity))) / purity

  return(n)

}

# Colour with saturation
cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
           "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
         rep("#FF4542", times = 100 - 8))

# Name top
names(cols)[(9:100)+1] = 9:100

dist2d = function(data_point,segment_point1=c(0,1),segment_point2=c(1,0), absolute=F) {
  v1 = segment_point1 - segment_point2
  v2 = data_point - segment_point1
  m = cbind(v1,v2)
  if(absolute) {d = abs(det(m))/sqrt(sum(v1*v1))} else {
    
    d = det(m)/sqrt(sum(v1*v1))
    d = d*-1
    
  }
  return(d)
} 

tertile_split = function(x) {
  res = ifelse(is.na(x), 
               yes = NA, 
               no = ifelse(x > quantile(x, probs = 1/3, na.rm = T), 
                           yes = ifelse(x > quantile(x, probs = 2/3, na.rm = T), 
                                        yes = "3rd", 
                                        no = "2nd"), 
                           no = "1st"))
  return(res)
}