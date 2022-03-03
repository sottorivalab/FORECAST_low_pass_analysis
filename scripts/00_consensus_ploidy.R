library(copynumber)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Plotting settings
run_full_pipeline = TRUE

# Individual plotting option
individual_colour   = c("ploidy", "gain_loss", "no_colour")[1]
plot_chromosome     = TRUE

source("scripts/00_general_functions.R")
source("runASCATlp.R")

# What are the ranges for lower and higher ploidy searches?
lower_range = seq(1.5,2.5,by=0.1)
highr_range = seq(3.5,4.5,by=0.1)

# Intervals for purity search
int          = 0.001
# Lowest purity allowed
min_purity   = 0.1
# Maximum Log2ratio allowed when fitting
max_lrr      = Inf
# If there is no fit, what purity do we use?
nofit_purity = 0.4

# Where is the low pass data?
pipeline_root   = "~/remotes/forecast_low_pass/"
# Add on a tag to the results for version testing
version_add_on = ""

# Create directories needed if don't exist
dir.create(paste0("results/0_min_het_best_fit_ploidy_search/calls",version_add_on), showWarnings = F)
dir.create(paste0("results/0_min_het_best_fit_ploidy_search/case_ploidy",version_add_on), showWarnings = F)
dir.create(paste0("results/0_min_het_best_fit_ploidy_search/plot_separate",version_add_on), showWarnings = F)
dir.create(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on), showWarnings = F)

# Get patient folders
patients = list.files(paste0(pipeline_root,"data"))
patients = patients[!patients %in% "README"]

# Read in the metrics analysis
per_sample_mets = read.table("results/2_metrics/metric_assessment_table.txt", 
                             header = T, stringsAsFactors = F)

# Make a blacklist of impure (by eye) or failed samples samples because they will bias the search
normal_samples = per_sample_mets$Sample[per_sample_mets$Assessment=="Failed" | per_sample_mets$Assessment=="Normal"]

# Collect the results
per_patient_solution = list()
med_dev_zeros = NULL

# Do it on a per patient basis
for(patient in patients) {
  
  # Find files via the Snakemake pipeline infratstructure
  files = list.files(paste0(pipeline_root,"data/",patient,"/QDNASEQ"), 
                     pattern = "_cna_segments.txt", full.names = T, recursive = T)

  # Read in data
  segs = lapply(files, function(f) read.table(f)[,1])
  
  # Find files
  files = list.files(paste0(pipeline_root,"data/",patient,"/QDNASEQ"), 
                     pattern = "_bins.txt", full.names = T, recursive = T)
  
  # Read in data
  bins = lapply(files, function(f) {
    
    dat = read.table(f, stringsAsFactors = F)
    
    dat$chromosome = factor(dat$chromosome, levels = c(1:22,"X","Y"))
    
    return(dat)
    
  })
  
  # For now we remove the sex chromosomes!
  autosome_index = bins[[1]]$chromosome %in% 1:22
  segs_auto = lapply(segs, function(i) i[autosome_index])
  
  # Bins of the autosomes
  bins_auto = lapply(bins, function(i) {
    
    res = i[autosome_index,]
    
    res$chromosome = as.numeric(res$chromosome)
    
    return(res)
    
  })
  
  # Read in the blacklist bins
  blacklist = readRDS("refs/Bins_to_blacklist_from_recurrence_analysis_20210723.rds")

  # Good bin index
  good_index = !bins_auto[[1]]$feature %in% paste0(blacklist$chr,":",blacklist$start,"-",blacklist$end)

  # Make test segs
  segs_test = lapply(segs_auto, function(i) i[good_index])

  # Extract names
  sample    = gsub("_dups_valid", "", unlist(lapply(bins_auto, function(b) colnames(b)[5])))
  sample_df = data.frame(do.call(rbind, strsplit(sample, split = "_")), stringsAsFactors = F)
  colnames(sample_df) = c("Patient", "Region", "DW")
  
  # Median deviation from zero
  med_dev_zero = sapply(1:length(segs_test), function(i) {
    
    if(sample[i] %in% per_sample_mets$Sample[per_sample_mets$Assessment=="Failed"]) {
      
      res = 0
      
    } else {res = median(abs(segs_test[[i]] - 0)^2)}
    
    return(res)
    
  })
  
  # Collect the maximum observed value in the patient
  med_dev_zeros = c(med_dev_zeros, max(med_dev_zero))
  
  # Essentially firstly call ploidy based on this measure
  high_ploidy = any(med_dev_zero > 0.003)
  
  # Info the user what the result is
  print(paste0(patient," median deviation=",max(med_dev_zero)))
  
  # Does it pass the threshold?
  if(high_ploidy) {ploid_search = highr_range} else {ploid_search = lower_range}
  
  # For each ploidy test fits
  ploidy_mean_dist = lapply(ploid_search, function(pld) {
    
    output = list()
    
    # Remove normal and failed samples from consideration
    seg_test = segs_test[!sample %in% normal_samples]
    
    # Default to diploid if there is nothing to tumour like
    if(length(seg_test)==0) {
      
      if(pld == 2) {
        
        output$int_dists = 0
        output$sol_dists = 0
        
      } else {
      
        output$int_dists = 100
        output$sol_dists = 100
        
        return(output)
      
      }
      
    } else {
    
      # What's the fit for this sample?
      ascat = lapply(seg_test, runASCATlp, fix_ploidy=pld, interval = int,
                     min_purity=min_purity, max_lrr = max_lrr, no_fit_psit = pld)
      
      # Normalise it by ploidy
      mat = do.call(rbind, lapply(ascat, function(s) s$CN)) / pld
      
      # Remove extreme LRR values from the dist measure across samples
      mat = mat[,!apply(do.call(cbind, seg_test), 1, max) > 1]
      
      # Calculate euclidean distances
      dists = as.vector(dist(mat))
      
      # Collect as output
      output$int_dists = dists
      output$sol_dists = unlist(lapply(ascat, function(i) {
        
        within_dist = i$allsol$dist[1]
        
        # If not fit, you give the measure of fit for a pure sample
        if(is.na(within_dist)) {within_dist = i$fit[rownames(i$fit)==1,]}
        
        return(within_dist)
        
      }))
      
      return(output)
      
    }
    
  })
  
  # Calculate means across tested ploidies
  inter_sample_dist = unlist(lapply(ploidy_mean_dist, function(i) mean(i$int_dists)))
  CN_solution_dist = unlist(lapply(ploidy_mean_dist, function(i) mean(i$sol_dists, na.rm=T)))

  # Collect as a dataframe
  df = data.frame(inter_sample = inter_sample_dist, in_sample = CN_solution_dist)
  
  # Scale it between zero and 1
  df_scaled = data.frame(apply(df, 2, range01))

  # What's the distance from 0,0 for each ploidy tested
  dist_from_zero = apply(df_scaled, 1, function(x) dist(rbind(x, rep(0, times = 2))))

  # Take the minimum
  ploidy_sol_index = which.min(dist_from_zero)
  
  # Now for making a plot to view the results
  plt.df = df_scaled
  plt.df$Solution = FALSE
  plt.df$Solution[ploidy_sol_index] = TRUE
  
  # Where is the solution?
  cols = c("FALSE" = "black", "TRUE" = "red")
  
  # Plot out the final analysis
  plt = ggplot(plt.df, aes(x = in_sample, y = inter_sample, col = Solution)) + geom_point() + 
    annotate("text", x = plt.df$in_sample, y = plt.df$inter_sample+0.03, label = ploid_search) + 
    scale_colour_manual(values = cols) +
    xlab("Within sample mean rel. distance") + ylab("Between samples rel. distance") +
    theme_bw()
  
  # Save it
  if(run_full_pipeline) {
    
    pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/",patient,"_02_scaled_distances.pdf"), width = 9)
    print(plt)
    dev.off()
    
  }
  
  # What's the answer?
  searched_ploidy = ploid_search[ploidy_sol_index]

  # Plot the ploidy v intersample fit
  if(run_full_pipeline) {
    
    pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/",patient,"_01_mean_intersample_distance_ploidy.pdf"), 
        width = 9)
    plot(ploid_search, unlist(lapply(ploidy_mean_dist, function(i) mean(i$int_dists))), type = "l",
         xlab = "Ploidy", ylab = "Mean Intersample Distance")
    dev.off()
    
  }
  
  # With this answer generate the final calls
  ascat = lapply(segs_auto, runASCATlp, fix_ploidy=searched_ploidy, interval = int,
                 min_purity=min_purity, max_lrr = max_lrr, no_fit_psit = searched_ploidy)
  names(ascat) = sample
  
  # Plot the fit of the segments with the solutions
  if(run_full_pipeline) {
    
    pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/",patient,"_05_fits.pdf"))
    
    for(s in 1:length(ascat)) {
      
      mean_dist = mean(ascat[[s]]$CN - ascat[[s]]$contCN)
      hist(ascat[[s]]$CN - ascat[[s]]$contCN, breaks = 100, xlim = c(-1,1), main = names(ascat)[s], xlab = "Distance")
      abline(v = mean_dist, lty = "dotted")
      
    }
    
    dev.off()
    
  }
  
  # If there is no solution make a fit with the nofit purity
  for(s in 1:length(ascat)) {
    
    if(nrow(ascat[[s]]$allsol)==0) {
    
      print(names(ascat)[s])
      
      lrrs  = segs_auto[[s]]
      gamma = 1
      
      # Get continuous copy number values for best fit
      rho = nofit_purity
      psi = (2*(1 - rho)) + (rho*searched_ploidy)
      n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
      
      # Make them integers
      n_int = round(n)
      
      # If we get minus states (i.e. small deletions in impure tumours)
      n_int[n_int<0] = 0
      
      ascat[[s]]$Purity = nofit_purity
      ascat[[s]]$Psi    = psi
      ascat[[s]]$contCN = n
      ascat[[s]]$CN     = n_int
    
    }}
  
  # Add on the X and Y chromosomes in accordance with the purity, ploidy fits
  for(s in 1:length(ascat)) {
    
    sex_segs = segs[[s]][bins[[s]]$chromosome %in% c("X", "Y")]
    
    n = callXchromsome(sex_lrrs = sex_segs, psi = ascat[[s]]$Psi, psit = ascat[[s]]$PsiT, purity = ascat[[s]]$Purity)
    
    ascat[[s]]$CN = c(ascat[[s]]$CN, round(n))
    ascat[[s]]$contCN = c(ascat[[s]]$contCN, n)
    
  }
  
  # Collect as a matrix
  mat = do.call(rbind, lapply(ascat, function(s) s$CN))
  
  # Colour with saturation
  cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
             "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
           rep("#FF4542", times = 100 - 8))
  
  # Name top
  names(cols)[(9:100)+1] = 9:100
  
  # Save the object for calling
  saveRDS(ascat, file = paste0("results/0_min_het_best_fit_ploidy_search/calls",version_add_on,"/",patient,"_ascat_objects.rds"))
  
  # Make a pdf of these calls together
  if(run_full_pipeline) {pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/",patient,"_03_new_calls.pdf"), width = 14)}
  
  plts = lapply(1:length(ascat), function(s) {
    
    sample_ascat = ascat[[s]]
    
    # Make a plotting dataframe
    plt.df = data.frame(genome.bin = 1:nrow(bins[[s]]), 
                        chr = bins[[s]]$chromosome, 
                        Log2ratio = bins[[s]][,5],
                        Call = as.character(sample_ascat$CN),
                        mean_segment = segs[[s]],
                        segment_col = "green")
    
    # Make plot
    p = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_hline(yintercept = c(-2,-1,0,1,2), lty = c("solid"), lwd = 0.2) +
      geom_point() +
      scale_colour_manual(values = cols) +
      scale_x_continuous(name = "Chromosomes", labels = c(1:22,"X","Y"), 
                         breaks = as.vector(c(1, cumsum(table(bins[[s]]$chromosome))[-24]) + 
                                              (table(bins[[s]]$chromosome) / 2))) + 
      geom_vline(xintercept = c(1, cumsum(table(bins[[s]]$chromosome))), lty = "dotted") +
      ggtitle(paste0("Low pass calls - ",names(ascat[s])," purity=",
                     sample_ascat$Purity,", psit = ",sample_ascat$PsiT,
                     ", ploidy = ",signif(mean(sample_ascat$CN[autosome_index]), digits = 3),")")) + 
      scale_y_continuous(limits=c(-2,2), oob=scales::squish) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18)) +
      geom_point(aes(y = mean_segment), color="#000000")
    
    print(p)
    
    return(p)
    
  })
  
  if(run_full_pipeline) {dev.off()}
  
  names(plts)   = names(ascat)
  rownames(mat) = names(ascat)
  
  # Also produce plots as individual files
  lapply(1:length(ascat), function(s) {
    
    sample_ascat = ascat[[s]]
    
    # Convert to losses and gains 
    ploidy = round(searched_ploidy)
    
    # Get em
    calls = sample_ascat$CN
    
    # It doesn't matter because it's a name
    cols = c(cols, "Neutral" = "grey", "Lost" = "#75bbfd", "Gain" = "#ff474c")
    
    if(individual_colour=="gain_loss") {
      
      # Baseline
      calls[which(sample_ascat$CN==ploidy)] = "Neutral"
      calls[which(sample_ascat$CN<ploidy)]  = "Lost"
      calls[which(sample_ascat$CN>ploidy)]  = "Gain"
      
      cols = c(cols, "Neutral" = "grey", "Lost" = "#75bbfd", "Gain" = "#ff474c")
      
    }
    
    if(individual_colour=="no_colour") {
      
      # Baseline
      calls = "Neutral"
      
    }
    
    # Make a plotting dataframe
    plt.df = data.frame(genome.bin = 1:nrow(bins[[s]]), 
                        chr = bins[[s]]$chromosome, 
                        Log2ratio = bins[[s]][,5],
                        Call = as.character(calls),
                        mean_segment = segs[[s]])
    
    # Max CN
    maxCN = 2
    minCN = -2
    
    # lines across
    lines_across = data.frame(x1 = 1,
                              y1 = minCN:maxCN,
                              x2 = nrow(bins[[s]]),
                              y2 = minCN:maxCN)
    
    # lines going up for chromosomes
    lines_vertical = data.frame(x1 = c(1, cumsum(table(bins[[s]]$chromosome))),
                                y1 = minCN,
                                x2 = c(1, cumsum(table(bins[[s]]$chromosome))),
                                y2 = maxCN)
    
    # Remove any chromosome labels due to congestion?
    chr_out = c(19,21)
    
    chrs_lab = c(1:22,"X","Y")
    chrs_lab[chr_out] = ""
    
    # Make plot
    p = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = lines_across, color="#00000080", lwd = 0.2) +
      scale_colour_manual(values = cols) + 
      scale_x_continuous(name = NULL, 
                         labels = chrs_lab,
                         breaks = as.vector(c(1, cumsum(table(bins[[s]]$chromosome))[-24]) + 
                                              (table(bins[[s]]$chromosome) / 2)),
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
    
    if(plot_chromosome) {
      pdf(paste0("results/0_min_het_best_fit_ploidy_search/plot_separate",version_add_on,"/",names(ascat)[s],"_higher_ploidy_calls.pdf"),
          height = 5, width = 24, useDingbats=FALSE)
      print(p)
      dev.off()
    }
    
    if(!plot_chromosome) {
      p = p + theme(axis.text.x = element_blank())
      ggsave(paste0("results/0_min_het_best_fit_ploidy_search/plot_separate",version_add_on,"/",names(ascat)[s],"_higher_ploidy_calls_no_chromo.png"), 
             plot = p,
             height = 5, width = 24)
    }
    
  })
  
  if(run_full_pipeline) {pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/",patient,"_04_heatmap.pdf"), width = 14)}
  
  # Let's also make a heatmap of the results
  hm = Heatmap(mat, cluster_rows = T, cluster_columns = F,
               col = colorRamp2(c(0, round(searched_ploidy), 8), c("blue", "white", "red")), 
               name = "CN", column_title = paste0(patient,"=",searched_ploidy))
  
  draw(hm)
  
  # Extract chromosome vector
  chrs = unlist(lapply(strsplit(bins[[1]]$feature, split = "[:]"), function(i) i[1]))
  
  #Chromosome separation positions
  chr.ends = cumsum(rle(chrs)$lengths)[-24]
  
  #Add lines
  for(boundary in c(0,1)) {
    
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
  for(boundary in chr.ends / ncol(mat)) {
    
    #Add the lines
    decorate_heatmap_body("CN", {
      grid.lines(c(boundary, boundary), c(0, 1), gp = gpar(lty = "dotted", lwd = 1))
    })
    
  }
  
  if(run_full_pipeline) {dev.off()}
  
  if(run_full_pipeline) {
    
    lapply(1:length(ascat), function(i) {
      
      CN_out = data.frame(ascat[[i]]$CN)
      
      colnames(CN_out) = names(ascat)[i]
      
      rownames(CN_out) = bins[[i]]$feature
      
      write.table(CN_out, file = paste0("results/0_min_het_best_fit_ploidy_search/calls",version_add_on,"/",sample_df[i,"Patient"],
                                        "_",
                                        sample_df[i,"Region"],"_",sample_df[i,"DW"],"_cna_ploidy_search_calls.txt"),
                  quote = F)
      
      # Calculate ploidy vector
      ploidy_expected = rep(1, times = length(ascat[[i]]$CN))
      ploidy_expected[!autosome_index] = 0.5
      ploidy_expected = round(searched_ploidy)*ploidy_expected
      ploidy_expected = round(ploidy_expected)
      
      metrics = data.frame(Sample = names(ascat)[i],
                           Purity = ascat[[i]]$Purity, 
                           PsiT = ascat[[i]]$PsiT, 
                           Ploidy = signif(mean(ascat[[i]]$CN[autosome_index]), digits = 4),
                           PGA = signif(length(which(ascat[[i]]$CN!=ploidy_expected)) / length(ascat[[i]]$CN), digits = 4)
                           )
      
      write.table(metrics, file = paste0("results/0_min_het_best_fit_ploidy_search/calls",version_add_on,"/",sample_df[i,"Patient"],
                                        "_",
                                        sample_df[i,"Region"],"_",sample_df[i,"DW"],"_cna_ploidy_search_metrics.txt"),
                  quote = F, row.names = F, sep = "\t")
      
    })
    
  }
  
  
  per_patient_solution[[patient]] = searched_ploidy
  
}

if(run_full_pipeline) {
  
  # Cohort ploidy range
  pdf(paste0("results/0_min_het_best_fit_ploidy_search/plots",version_add_on,"/FORECAST_cohort_ploidy_range.pdf"))
  hist(unlist(per_patient_solution), breaks = 100, xlim = c(1.5,4.5),
       main = "PsiT solutions", xlab = "Tumour Ploidies", ylab = "Frequency")
  dev.off()
  
}

# Make a ploidy list per case
per_case_psit = data.frame(Case = names(per_patient_solution),
                           PsiT = unlist(per_patient_solution))
rownames(per_case_psit) = NULL

if(run_full_pipeline) {
  
  # Write out the psit 
  write.table(per_case_psit, row.names = F, col.names = T, quote = F, sep = "\t",
              file = paste0("results/0_min_het_best_fit_ploidy_search/case_ploidy",version_add_on,"/patient_psit_solutions.txt"))
  
}
