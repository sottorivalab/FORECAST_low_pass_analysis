# Script for detecting amplifications in a bespoke manner using log2ratio values 

# Library
library(ggplot2)
library("biomaRt")
library("copynumber")

# Functions
source("scripts/00_general_functions.R")

# Parameters
amp_min = 1
amp_max = 20
quant   = 0.99
minlrr  = 0.58
z_score = 3

# Setting
mode = c("quantile", "zscore", "hard")[2]

# Where is the low pass data?
pipeline_root = "~/remotes/forecast_low_pass/"

# Read in the chromosome arm boundaries
arms = read.table("refs/chrArmBoundaries_hg38.txt", header = T)

# Read in the prognostic genes from the Human Protein Atlas
prog_genes = read.table("refs/prognostic_prostate.tsv", header = T, stringsAsFactors = F, sep = "\t")

# Set up BiomaRt
listMarts(host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)

# Query_list
query = c("MDM2",
          "MYCN",
          "AR",
          "TERC", 
          "MYBL3", 
          "HRAS", 
          "PI3KCA", 
          "JUNB", 
          "LAMC2", 
          "RAF1", 
          "MYC", 
          "GARP", 
          "SAS", 
          "FGFR1", 
          "PGY1", 
          "MYCL1", 
          "MYB", 
          "FGR", 
          "HER2", 
          "CCND1", 
          "FOXA1",
          prog_genes$Gene)

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Find the log2ratio files
per_patient_samples = lapply(patients, FORECAST_list.files, pattern = "bins.txt")

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", header = T)

# Get the failed samples
failed_samples  = paste0(per_sample_mets[per_sample_mets$Assessment=="Failed","Patient"],"_",
                         per_sample_mets[per_sample_mets$Assessment=="Failed","Region"])

# Name the patient entries
names(per_patient_samples) = patients

# List for recurrence analysis
genes_cohort = NULL

# Collect for per patient summary
per_patient_summary = makeSummaryList(patients)

# Case name
for(case in patients) {
  
  # Patient files to read in 
  patient_files = per_patient_samples[[case]]
  
  # Remove failed samples
  patient_files = patient_files[!gsub("_bins.txt", "", basename(patient_files)) %in% failed_samples]
  
  # Read in the data as a list
  samples = FORECAST_read.files(patient_files)
  
  # Gather here the chromosome position
  chr_pos = samples[[1]][,2:4]
  
  # Process the samples
  lrrs = lapply(samples, function(l) {
    
    cn = gsub("_dups_valid", "", colnames(l)[5])
    
    l = cbind(l[,5])
    
    colnames(l) = cn
    
    return(l)
    
  })
  
  # Make a matrix out of it
  lrrs = do.call(cbind, lrrs)
  
  # Is it just a single sample?
  single_sample = ncol(lrrs)==1
  
  # Do single pcf if just one sample, if not do joint
  if(single_sample) {
    
    mssegs = pcf(data.frame(chr = chr_pos$chromosome, 
                            pos = chr_pos$start, 
                            lrrs), 
                 arms = getPQ(chr_pos, arms), 
                 gamma = 15, fast = FALSE)
    
  } else {
  
    # Segment the data using joint segmentation
    mssegs = multipcf(data.frame(chr = chr_pos$chromosome, 
                                 pos = chr_pos$start, 
                                 lrrs), 
                      arms = getPQ(chr_pos, arms), 
                      gamma = 15, fast = FALSE)
    
  }
  
  per_patient_summary[[case]] = makeSummaryList(gsub("_bins.txt", "", basename(patient_files)))
    
  for(sf in patient_files) {
    
    # Get base name of sample
    base_name = gsub("_bins.txt", "", basename(sf))
    
    # Read the bins data
    bins = read.table(sf, 
                      sep = "\t", header = T)
    
    # Dictate levels
    bins$chromosome = factor(bins$chromosome, levels = c(1:22, "X", "Y"))
    
    # Make log ratio column name generic
    colnames(bins)[5] = "log2ratio"
    
    # Make lrrs a separate object for manipulation
    lrrs = bins[,5]
    
    # Normalise for male sex chromosomes
    norm_lrrs = male_normalise_lrr(logr = lrrs, chr_vec = bins$chromosome)
    
    # Calculate upper threshold to use
    if(mode == "quantile") {
      amp_threshold = max(minlrr, quantile(norm_lrrs, probs = quant))
    }
    if(mode == "zscore") {
      amp_threshold = (z_score*sd(norm_lrrs)) + mean(norm_lrrs)
    }
    if(mode == "hard") {
      amp_threshold = minlrr
    }
    
    # Get column to subset
    if(single_sample) {sub_col = "mean"} else {sub_col = base_name}
    
    # Male normalise segment output
    seg_male_norm = male_normalise_lrr(logr = mssegs[,sub_col], mssegs$chrom)
    
    # Get amplified bins from joint segmentation
    amp_df = mssegs[seg_male_norm > amp_threshold,c("chrom","start.pos","end.pos","n.probes")]
    
    # Original ends to do a check
    orig_end = amp_df$end.pos
    
    # Make end position the actually end of the bin location
    amp_df$end.pos = chr_pos$end[match(paste0(amp_df$chrom,"_",amp_df$start.pos), paste0(chr_pos$chromosome,"_",chr_pos$start)) + (amp_df$n.probes - 1)]
    
    # Throw error if not all correct
    if(!all((amp_df$end.pos - orig_end) < 500000)) {stop("Not all bin ends are correct!!!")}
    
    # Subset amplification candidates
    final_amps = amp_df[amp_df$n.probes>=amp_min & amp_df$n.probes <= amp_max,]
    
    # Get amplification coordinates
    amp_coords = lapply(1:nrow(final_amps), function(r) {
      
      row = final_amps[r,]
      
      c = row[,"chrom"]
      s = row[,"start.pos"]
      e = row[,"end.pos"]
      
      #loc = data.frame(chr = chr, start = str, end = end)
      loc = paste0(c,":",s,":",e)
      
    })
    
    Sys.sleep(1)
    # Get overlapping genes from biomaRt
    results=getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                  filters = c("chromosomal_region", "biotype"),
                  values = list(chromosomal_region=amp_coords, biotype="protein_coding"), 
                  mart = ensembl)
    Sys.sleep(1)
    
    # Add to collection to assess subclonality
    per_patient_summary[[case]][[base_name]] = results
    
    # Write out per case gene overlaps
    write.table(results, 
                file = paste0("results/3_amplification_detection/gene_overlaps/",base_name,"_genes_overlapping_with_amps.csv"),
                sep = ",",
                quote = F,
                row.names = F)
    
    # Add a column for the amplifications
    bins$Amplification = FALSE
    
    # Calculate indices which are amplified
    if(nrow(final_amps)==0) {
      
      amp_index = NA
      
    } else {
      
      amp_index = unlist(apply(final_amps, 1, function(r) {
        
        ac = r[1]
        as = as.numeric(r[2])
        ae = as.numeric(r[3])

        s = which(bins$chromosome==ac & bins$start==as)
        e = which(bins$chromosome==ac & bins$end==ae)
        
        ind = s:e
        
        return(ind)
        
      }))
    
    }
    
    # Is it considered amplified?
    bins$Amplification[amp_index] = TRUE
    
    bins$index = 1:nrow(bins)
    
    mb_min = amp_min*500 / 1000
    mb_max = amp_max*500 / 1000
    
    # Overlapping genes
    genes_overlapping = results$hgnc_symbol[results$hgnc_symbol %in% query]
    
    amp_genes = "-"
    
    if(length(genes_overlapping)>0) {
      
      amp_genes = paste0(paste0(genes_overlapping,"+"), collapse = ", ")
      
      # Re-search position for annotation
      key_gene_pos = getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                           filters = c("hgnc_symbol"),
                           values = list(hgnc_symbol=genes_overlapping), 
                           mart = ensembl)
    
      # Add middle of gene
      key_gene_pos$mid_position = ((key_gene_pos$end_position - key_gene_pos$start_position)/2) + key_gene_pos$start_position
      
      # Find the overlapping bins
      key_gene_pos$bindex = unlist(lapply(1:nrow(key_gene_pos), function(r) {
        
        bin = which(bins$chromosome==key_gene_pos$chromosome_name[r] & 
                      bins$start <= key_gene_pos$mid_position[r] &
                      bins$end > key_gene_pos$mid_position[r])
        
        # If the middle doesn't overlap try the start
        if(length(bin)==0) {
          
          bin = which(bins$chromosome==key_gene_pos$chromosome_name[r] & 
                        bins$start <= key_gene_pos$start_position[r] &
                        bins$end > key_gene_pos$start_position[r])
          
        }
        
        # If the start doesn't overlap try the end
        if(length(bin)==0) {
          
          bin = which(bins$chromosome==key_gene_pos$chromosome_name[r] & 
                        bins$start <= key_gene_pos$end_position[r] &
                        bins$end > key_gene_pos$end_position[r])
          
        }
        
        # OK if neither work find the nearest neighbour to the left
        if(length(bin)==0) {
          
          left  = max(which(bins$chromosome==key_gene_pos$chromosome_name[r] & bins$end <= key_gene_pos$mid_position[r]))
          right = left+1
          
          bin = c(left,right)[which.min(c(abs(bins[left,"end"] - key_gene_pos$mid_position[r]), abs(bins[right,"start"] - key_gene_pos$mid_position[r])))]
          
        }
        
        return(bin)
        
      }))
      
      # Create a random height
      key_gene_pos$height = seq(1.825,2.475, length.out = nrow(key_gene_pos))
    
    }
    
    # Collect for cohort
    genes_cohort = c(genes_cohort, results$hgnc_symbol)
    
    # Make plot
    p = ggplot(bins, aes(x = index, y = log2ratio, col = Amplification)) +
      geom_hline(yintercept = c(-2,-1,0,1,2,3), lty = c("solid"), lwd = 0.2) +
      geom_point() +
      scale_colour_manual(values = c(c("FALSE" = "grey", "TRUE" = "#E69F00"))) +
      scale_x_continuous(name = "Chromosomes", labels = c(1:22,"X","Y"), 
                         breaks = as.vector(c(1, cumsum(table(bins$chromosome))[-24]) + (table(bins$chromosome) / 2))) + 
      geom_vline(xintercept = c(1, cumsum(table(bins$chromosome))), lty = "dotted") +
      scale_y_continuous(limits=c(-2,3), oob=scales::squish) + 
      ggtitle(paste0("Low pass amplication calls - ",base_name,
                     " (threshold=",round(amp_threshold, digits = 2),
                     ", range=",mb_min,"-",mb_max,"Mb",
                     ", ",amp_genes,")")) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18))
    
    if(length(genes_overlapping)>0) {
      
      p = p + annotate("segment", 
                       x = key_gene_pos$bindex-200, 
                       xend = key_gene_pos$bindex, 
                       y = key_gene_pos$height, yend = bins[key_gene_pos$bindex,"log2ratio"], 
                       colour = "black", size=0.5, alpha=1) +
        annotate("text", x = key_gene_pos$bindex-360, y = key_gene_pos$height, 
                 label = paste0("italic(",key_gene_pos$hgnc_symbol,")"), parse = TRUE)
      
    }
    
    p
    
    ggsave(paste0("results/3_amplification_detection/amp_plots/",base_name,"_amplification_detection.pdf"), height = 5, width = 12)
    
  }
  
}

# Write out the genes that appear more than once
write.table(data.frame(table(genes_cohort)), file = "results/3_amplification_detection/gene_overlaps/0_cohort_summary_gene_frequencies_analysis.csv",
            sep = ",",
            quote = F,
            row.names = F)

# Summarise subclonality
subclonality_summary = lapply(per_patient_summary, function(p) {
  
  collapsed = do.call(rbind, p)
  collapsed = collapsed[collapsed$hgnc_symbol!="",]
  res = table(collapsed$hgnc_symbol) / length(p)
  
  return(res)
  
})

# Save subclonality assessment
saveRDS(subclonality_summary, "results/3_amplification_detection/subclonality_assessment.rds")
