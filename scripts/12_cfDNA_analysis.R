library(copynumber)
library(ggplot2)
library(cowplot)
source("scripts/00_general_functions.R")
source("runASCATlp.R")

# These are the arms of hg38
arms       = read.table("chrArmBoundaries_hg38.txt", header = T)

# Simplistic copy number calling
cfDNA_patients = c("")

########################################################################################################
########### Perform multisample segmentation and call copy number alterations using lpASCAT ############
########################################################################################################

# Find the cfDNA samples
files = lapply(cfDNA_patients, function(patient) {
  
  dirs = list.files(paste0("~/remotes/forecast_low_pass/data/",patient,"/insert_size_subset"), 
                    full.names = T, recursive = T, pattern = "_bins.txt")
  
  dirs[grep("cfDNA", dirs)]
  
})

# Read the files in and store as list
per_patient_lrrs = lapply(files, function(fi) {
  
  lapply(fi, function(f) read.table(f, header = T, stringsAsFactors = F))
  
})

# Get the index of bins
chr_pos = per_patient_lrrs[[1]][[1]][,2:4]
bin_row = per_patient_lrrs[[1]][[1]][,1]

# chr_pos
chr_pos$chromosome = factor(chr_pos$chromosome, levels = c(1:22,"X","Y"))

# Take only the lrr data
per_patient_lrrs = lapply(per_patient_lrrs, function(p) {
  
  lapply(p, function(l) {
  
    cn = gsub("_dups_valid", "", colnames(l)[5])
    
    l = cbind(l[,5])
    
    colnames(l) = cn
    
    return(l)
  
  })
  
})

# Make it a matrix
per_patient_lrr_mat = lapply(per_patient_lrrs, function(p) do.call(cbind, p))

# Perform multi-sample segmentation calling
per_patient_ms_segs = lapply(per_patient_lrr_mat, function(p_mat) {
  
  if(ncol(p_mat) > 1) {

    # Perform the multiregion segmentation
    per_patient_multisample_seg = multipcf(data.frame(chr = chr_pos$chromosome, 
                                                      pos = chr_pos$start, 
                                                      p_mat), 
                                           arms = getPQ(chr_pos, arms), 
                                           gamma = 10, fast = FALSE)
    
    
    
  } else {
    
    per_patient_multisample_seg = pcf(data.frame(chr = chr_pos$chromosome, 
                                                 pos = chr_pos$start, 
                                                 p_mat), 
                                      arms = getPQ(chr_pos, arms), 
                                      gamma = 10, fast = FALSE)
    
    colnames(per_patient_multisample_seg)[7] = per_patient_multisample_seg$sampleID[1]
    per_patient_multisample_seg$sampleID = NULL
  
  }
  
  # Expand it out per sample
  per_patient_ms_segs = lapply(1:(ncol(per_patient_multisample_seg) - 5), function(s) {
    
    rep(per_patient_multisample_seg[,s+5], times = per_patient_multisample_seg[,5])
    
  })
  
  # Name the elements in the list
  names(per_patient_ms_segs) = colnames(per_patient_multisample_seg)[6:ncol(per_patient_multisample_seg)]
  
  return(per_patient_ms_segs)
  
})

# Run lpASCAT
per_patient_ms_ascat = lapply(1:length(per_patient_ms_segs), function(i) {
  
  p = per_patient_ms_segs[[i]]
  
  lapply(1:length(p), function(s) {
    
    # Write out the segments for use in the manual ploidy choosing script
    samples = names(p)
    
    mid_pld = 3.1
    expand  = 1.6
    mp      = 1
    
    if(any(grepl("XXXXX", samples))) {
      
      mid_pld = 4.35
      expand  = 0.35
      
    }
    
    autosome_index = chr_pos$chromosome %in% 1:22
    
    bins_auto = per_patient_lrr_mat[[i]][autosome_index,s]
    segs_auto = p[[s]][autosome_index]
    
    sn = samples[s]
    
    ps = F
    pr = 1000
    pp = 1000
    
    if(sn == "XXXXX") {
      
      ps = T
      pr = 0.07
      pp = 4.41
      
    }
    
    res = runASCATlp(lrrs = segs_auto, fix_ploidy = mid_pld, pad_ploidy = expand,
                     interval = 0.01, min_purity=0.01, max_lrr = Inf, no_fit_psit = 2, preset = ps, 
                     preset_purity = pr, preset_ploidy = pp, max_purity = mp)
    
    sex_segs = p[[s]][!autosome_index]
    
    n = callXchromsome(sex_lrrs = sex_segs, psi = res$Psi, psit = res$PsiT, purity = res$Purity)
    
    res$CN = c(res$CN, ifelse(round(n) < 0, 0, round(n)))
    res$contCN = c(res$contCN, n)
    
    res$segs       = p[[s]]
    res$bins       = per_patient_lrr_mat[[i]][,s]
    res$samplename = sn
    
    return(res)
    
  })
  
})

# Make plots per sample
lapply(per_patient_ms_ascat, function(i) {
  
  lapply(i, function(j) {
    
    cna_data = j
    
    sample_name = cna_data$samplename
    
    cn_output   = data.frame(cna_data$CN)
    colnames(cn_output) = sample_name
    
    # Make a plot dataframe
    plt.df = data.frame(genome.bin = 1:length(cna_data$bins),
                        Chromosome = chr_pos$chromosome,
                        Log2ratio = cna_data$bins,
                        mean_segment = cna_data$segs,
                        Call = as.factor(cna_data$CN))
    
    # Max CN
    maxCN = 2
    minCN = -2
    
    # lines across
    lines_across = data.frame(x1 = 1,
                              y1 = minCN:maxCN,
                              x2 = length(cna_data$bins),
                              y2 = minCN:maxCN)
    
    # lines going up for chromosomes
    lines_vertical = data.frame(x1 = c(1, cumsum(table(plt.df$Chromosome))),
                                y1 = minCN,
                                x2 = c(1, cumsum(table(plt.df$Chromosome))),
                                y2 = maxCN)
    
    # Remove any chromosome labels due to congestion?
    chr_out = c(19,21)
    
    chrs_lab = c(1:22,"X","Y")
    chrs_lab[chr_out] = ""
    
    # Make the plot
    p = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_hline(yintercept = c(-2,-1,0,1,2), lty = c("solid"), lwd = 0.2) +
      geom_point() +
      scale_colour_manual(values = cols) +
      scale_x_continuous(name = "Chromosomes", labels = c(1:22,"X","Y"), 
                         breaks = as.vector(c(1, cumsum(table(plt.df$Chromosome))[-24]) + 
                                              (table(plt.df$Chromosome) / 2))) + 
      geom_vline(xintercept = c(1, cumsum(table(plt.df$Chromosome))), lty = "dotted") +
      ggtitle(paste0("Low pass calls - ",sample_name,", purity=",
                     cna_data$Purity,", psit = ",cna_data$PsiT)) + 
      scale_y_continuous(limits=c(-2,2), oob=scales::squish) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18)) +
      geom_point(aes(y = mean_segment), color="#000000")

    ggsave(filename = paste0("results/10_cfDNA/",sample_name,"_multiregion_seg_cna_profile.png"),
           plot = p, width = 15, height = 5)
    
    write.table(cn_output, file = paste0("results/10_cfDNA/",sample_name,"_multiregion_seg_cna_calls.txt"),
                row.names = T, col.names = T, quote = F)
    
    # Calculate ploidy vector
    autosome_index = chr_pos$chromosome %in% 1:22
    
    ploidy_expected = rep(1, times = length(cna_data$CN))
    ploidy_expected[!autosome_index] = 0.5
    ploidy_expected = median(cna_data$CN)*ploidy_expected
    ploidy_expected = round(ploidy_expected)
    
    # Get output metrics
    metrics = data.frame(Sample = sample_name, Purity = cna_data$Purity, PsiT = cna_data$PsiT, 
                         Ploidy = signif(mean(cna_data$CN[autosome_index]), digits = 4), 
                         PGA = signif(length(which(cna_data$CN!=ploidy_expected)) / length(cna_data$CN), digits = 4))
    
    if(median(cna_data$CN) != round(metrics$Ploidy)) {warnings("Median CN and mean ploidy mismatch!")}
    
    write.table(metrics, file = paste0("results/10_cfDNA/",sample_name,"_multiregion_seg_metrics.txt"),
                row.names = F, col.names = T, quote = F, sep = "\t")
    
    # Make plot
    p2 = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_across, color="#00000080", lwd = 0.2) +
      scale_colour_manual(values = cols) + 
      scale_x_continuous(name = NULL, 
                         labels = chrs_lab,
                         breaks = as.vector(c(1, cumsum(table(plt.df$Chromosome))[-24]) + 
                                              (table(plt.df$Chromosome) / 2)),
                         expand = c(0.01,0)) +
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_vertical, color="#0000001A", lwd = 0.2) +
      scale_y_continuous(limits=c(-2,2), oob=scales::squish, expand = c(0,0)) + 
      geom_point(aes(size = 2)) +
      geom_point(aes(y = mean_segment), color="#000000") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.x = element_text(size = 18),
            #axis.text.x = element_blank(),
            axis.text.x = element_text(size = 40, angle = 45, vjust = 1, hjust=1),
            axis.title.y = element_blank(),
            axis.text.y = element_text(size = 60),
            legend.position = "none",
            plot.margin = margin(0.4, 0.1, 0.4, 0.1, "in"))
    
    pdf(paste0("results/10_cfDNA/",sample_name,"_multiregion_seg_cna_profile_separate.pdf"),
        height = 5, width = 24, useDingbats=FALSE)
    print(p2)
    dev.off()
    
  })
  
})

########################################################################################################
########### Create input for MEDICC2 analysis jointly with primary samples #############################
########################################################################################################

library(reshape2)

# Functions
source("scripts/00_general_functions.R")

# cfDNA patients
cfDNA_patients = c("")

medicc_root = "~/remotes/medicc_trees/"

min_pga = 0.01
min_bin_num  = 0

# Get patient folders
patients = cfDNA_patients

# Get all files
call_files = list.files("results/0_min_het_best_fit_ploidy_search/calls/", 
                        pattern = "_cna_ploidy_search_calls.txt",
                        full.names = T)

call_files = c(call_files, list.files("results/10_cfDNA", pattern = "IS90150_multiregion_seg_cna_calls.txt", full.names = T))

# Remove the impure cfDNA sample
call_files = call_files[!grepl("XXXXX", call_files)]

# Find the call files
per_patient_samples = lapply(patients, function(p) {call_files[grep(p, x = call_files)]})

# Read in the data
per_patient_calls = lapply(per_patient_samples, FORECAST_read.files)

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", header = T, stringsAsFactors = F)

# Make a blacklist of impure samples
flat_profiles = per_sample_mets$PGA < min_pga | 
  per_sample_mets$Assessment=="Failed" | 
  grepl("_N[0-9]{1,2}_", per_sample_mets$Sample) | 
  per_sample_mets$Sample==""
flat_profiles = per_sample_mets$Sample[flat_profiles]

# Get indices of the tumour samples
per_patient_calls = lapply(per_patient_calls, function(p) {
  
  t_ind = which(unlist(lapply(p, function(n) {!removeBamNameEnding(colnames(n)) %in% flat_profiles})))
  
  return(p[t_ind])
  
})

per_patient_mats        = lapply(per_patient_calls, function(p) do.call(cbind, p))
names(per_patient_mats) = patients

n_samples = unlist(lapply(per_patient_mats, ncol))

for(p in patients) {
  
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
    
    print(tail(chr_pos))
    
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
                     size = run_len$lengths,
                     m[cumulative_len,])
    
    # Filter for size
    res = res[res$size >= min_bin_num,]
    res$size = NULL
    
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
    write.table(res, file = paste0(medicc_root,"forecast_cfDNA/",p,"_regions.tsv"),
                sep = "\t", row.names = F, quote = F)
    
  }
  
}

########################################################################################################
########### Create plots from MEDICC2 analysis with cfDNA and primary samples ##########################
########################################################################################################

library(ggtree)
library(circlize)

# cfDNA patients
cfDNA_patients = c("")

lapply(cfDNA_patients, function(patient) {

  # Read in the file
  tree = read.tree(paste0(medicc_root,"forecast_cfDNA/",patient,"_regions_final_tree.new"))
  
  # Nick's colour scale for gleason
  gleason_cols = colorRamp2(c(1,2,3,4,5), colors = c("#009E73", "#F0E442", "#D55E00", "#5D3A9B", "#000000"))
  
  # Get per region gleason scores
  gleason_cont = read.csv("results/8_collect_metrics/Gleason_per_sample_digital_pathology.csv")
  
  # Get patient gleason
  p_gleason = gleason_cont[gleason_cont$Patient==patient,]
  
  # Make a short df of tree samples
  sample_gleason = data.frame(Sample = tree$tip.label)
  
  old = sample_gleason
  
  # Get cont_gleason
  sample_gleason = merge(sample_gleason, p_gleason[,c("Sample", "Grade.Group", "Cont_Gleason")], 
                         all.x = T, sort = F)
  
  sample_gleason = sample_gleason[match(old$Sample, sample_gleason$Sample),]
  
  print(all(old$Sample == sample_gleason$Sample))
  
  # Get the colour
  sample_gleason$colour = gleason_cols(sample_gleason$Grade.Group)
  
  # Make NAs some colour
  sample_gleason$colour[is.na(sample_gleason$colour)] = "grey"
  sample_gleason$colour[sample_gleason$Sample=="diploid"] = "#009E73"
  sample_gleason$colour[grepl("cfDNA", sample_gleason$Sample)] = "#ff000d"
  
  # Replace the patient and DW1
  tree$tip.label = gsub("_", " ", gsub("_DW1", "", gsub(paste0(patient,"_"), "", tree$tip.label)))
  tree$tip.label = gsub(" IS90150", "", tree$tip.label)
  
  # Plot
  p = ggtree(tree, size = 1) + geom_tiplab(hjust = -0.5, size = 8) + 
    geom_tippoint(size = 6, col = sample_gleason$colour) +
    ggtitle(patient) + theme(title=element_text(size=20, colour = "gray20"), 
                                                 plot.title = element_text(hjust = 0.5))
  
  # Add events labels?
  p = p + geom_label(data = p$data[p$data$branch.length > 0,], aes(x = branch, label = branch.length),
                     label.size = 0, nudge_y = 0.2)
  
  # Customise width?
  p2 = p
  p  = p2 + ggplot2::xlim(0, max(ape::node.depth.edgelength(tree)) * 1.3)
  ggsave(plot = p, filename = paste0("results/10_cfDNA/",patient,"_medicc2_cfDNA_trees.png"), 
         height = 6, width = 8)

})

########################################################################################################
########### Compare metrics and PGA specifically to the primary samples in the cohort ##################
########################################################################################################

# What are the cfDNA metric files
metric_files = list.files("results/10_cfDNA/", pattern = "_IS90150_multiregion_seg_metrics.txt", full.names = T)

# Read them in
metrics_df   = do.call(rbind, lapply(metric_files, read.table, header = T, stringsAsFactors = F))

# Read in all that lovely data
collected = read.table(file = "results/8_collect_metrics/FORECAST_genomic_metrics.txt", 
                       sep = "\t", stringsAsFactors = F, header = T)

# Subset now by the minimum number of samples required
collected = collected[collected$Samples >= 3,]

# Add a column for the cohort id
plot_df        = collected[,c("Patient", "mPGA")]
plot_df$cohort = "Primary"
colnames(plot_df) = c("Sample", "PGA", "cohort")

# Attach
metrics_df$cohort = "Relapse"
metrics_df = metrics_df[!metrics_df$Sample %in% c(""),]
plot_df = rbind(plot_df, metrics_df[,c("Sample", "PGA", "cohort")])

options(scipen = 10)
pvalue = signif(wilcox.test(plot_df$PGA[plot_df$cohort=="Relapse"], plot_df$PGA[plot_df$cohort=="Primary"])$p.value, 
                digits = 3)

# Plot the ctDNA PGAs
p = ggplot(plot_df, aes_string(x = "cohort", y = "PGA")) + geom_violin(lwd = 0.25) + 
  geom_jitter(width=0.1, lwd = 0.25) + 
  xlab("") + ylab("PGA") + ylim(0,1) + 
  geom_segment(aes(x = 0.992, y = 0.8, xend = 2.005, yend = 0.8), lwd = 0.25) +
  geom_segment(aes(x = 1, y = 0.8 - 0.01, xend = 1, yend = 0.8), lwd = 0.25) +
  geom_segment(aes(x = 2, y = 0.8 - 0.01, xend = 2, yend = 0.8), lwd = 0.25) +
  annotate("text", label = paste0("p=",pvalue), x = 1.5, y = 0.8 + 0.05, 
           size = 2) +
  annotate("text", label = paste0("n=",table(plot_df$cohort)), x = 1:2, y = 0.9, 
           size = 2) +
  theme_cowplot() + 
  theme(axis.text = element_text(size = 7), 
        axis.title = element_text(size = 7),
        axis.line = element_line(size = 0.25),
        axis.ticks = element_line(size = 0.25),
        axis.title.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), units = "mm"))
ggsave(plot = p, "results/10_cfDNA/PGA_comparison_cfDNA_DBs.pdf", 
       height = 47.33, width = 49.3, units = "mm", useDingbats = F)

# Now same for metrics df
metrics_df$Patient = unlist(lapply(strsplit(metrics_df$Sample, split = "_"), function(i) i[1]))

# Plot the ctDNA purities
p = ggplot(metrics_df, aes_string(x = 0, y = "Purity")) + geom_violin() + geom_jitter(width=0.1, aes(col = Patient)) +
  xlab("Relapse cfDNA") + ylab("Proportion of ctDNA in plasma") + scale_y_continuous(trans='log10') +
  theme_cowplot() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(plot = p, "results/10_cfDNA/Purity_cfDNA_distribution.png", height = 3, width = 10/3.5)

# Output info
write.csv(file = "results/10_cfDNA/cfDNA_analysis_metrics.csv", metrics_df, quote = F, row.names = F)
