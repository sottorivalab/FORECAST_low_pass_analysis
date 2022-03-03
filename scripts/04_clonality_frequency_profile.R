# Clonality frequency plot 
# Based on a plot originally created by Timon Heide (github.com/T-Heide/epicc_analysis_pub)

library(magrittr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(GenomicRanges)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
genome = BSgenome.Hsapiens.UCSC.hg38

# Functions
source("scripts/00_general_functions.R")

cutoff_clonal = 0.95 # fraction of samples in a case to be gained or lossed to consider site clonal
cutoff_subclonal = 0.0 # fraction of samples in a case to be gained or lossed to consider site clonal

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

# Read in the patient ploidies
patient_ploidies = read.table("results/0_min_het_best_fit_ploidy_search/case_ploidy/patient_psit_solutions.txt", 
                              header = T, stringsAsFactors = F)
patient_ploidies$Ploidy = round(patient_ploidies$PsiT)

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

# Make windows object from first sample
windows       = extractChrRows(rownames(per_patient_calls[[1]][[1]]))
windows$name  = paste0("win",1:nrow(windows))

# Reformat for plotting use
per_patient_calls = lapply(per_patient_calls, function(p) {
  
  # Perform reformating
  reformated = lapply(p, function(s) {
    
    # Reformating
    x = windows
    x$sample_barcode = removeBamNameEnding(colnames(s))
    x$value = s[,1]
    x$patient = extractPatientName(x$sample_barcode)
    
    return(x)
    
  })
  
  # Collect them as one
  reformated = do.call(rbind, reformated)
  
  return(reformated)
  
})

# Make one big dataframe
data_per_window = do.call(rbind, per_patient_calls)

# get cn gains and losses (relative to median cn):
id_vars = list(data_per_window$patient, data_per_window$name) # patient and window

# Make an autosome only data_per_window
data_per_window_autos = data_per_window[data_per_window$chr %in% c(1:22),]

median_cn_vector = patient_ploidies$Ploidy[match(data_per_window$patient, patient_ploidies$Case)]

median_cn_vector = ifelse(data_per_window$chr%in%c("X","Y"), 
                          median_cn_vector/2, 
                          median_cn_vector)

gained = data_per_window$value > median_cn_vector
clonal_status_gains = t(tapply(gained, id_vars, mean))

lossed = data_per_window$value < median_cn_vector
clonal_status_losses = t(tapply(lossed, id_vars, mean))

# censor all incomplete windows: 
stopifnot(all.equal(rownames(clonal_status_gains), rownames(clonal_status_losses)))
wh = apply(is.na(cbind(clonal_status_gains, clonal_status_losses)), 1, any)
clonal_status_gains[wh,] = NA
clonal_status_losses[wh,] = NA

################################################################################


data_per_window_patient_summary_gains = 
  cbind(windows, clonal_status_gains[windows$name,]) %>%
  reshape2::melt(id.vars=colnames(windows))


data_per_window_patient_summary_losses = 
  cbind(windows, clonal_status_losses[windows$name,]) %>%
  reshape2::melt(id.vars=colnames(windows)) %>% 
  mutate(value=-value)

chr_centers = 
  filter(D3GB::GRCh38.bands) %>%
  split(., .$chr) %>% 
  sapply(function(x) max(x$end)/2)

centromer_position = 
  filter(D3GB::GRCh38.bands, score=="acen") %>%
  split(., .$chr) %>% 
  sapply(function(x) mean(range(c(x$end,x$start))))

data_per_window_summary = 
  windows %>% 
  cbind(clonal_gains=apply(clonal_status_gains[windows$name,] > cutoff_clonal, 1, sum, na.rm=TRUE)) %>%
  cbind(subclonal_gains=apply(clonal_status_gains[windows$name,] > cutoff_subclonal, 1, sum, na.rm=TRUE)) %>%
  cbind(clonal_losses=apply(clonal_status_losses[windows$name,] > cutoff_clonal, 1, sum, na.rm=TRUE)) %>%
  cbind(subclonal_losses=apply(clonal_status_losses[windows$name,] > cutoff_subclonal, 1, sum, na.rm=TRUE)) 

ylm = ceiling(max(c(data_per_window_summary$subclonal_gains, data_per_window_summary$subclonal_losses)) / 10) * 10

plot_cna_recurence=
  data_per_window_summary %>%
  mutate(chr=gsub("chr","",chr)) %>%
  mutate(pos=(start+end)/2-chr_centers[as.character(chr)]) %>%
  mutate(chr=factor(chr, c(1:22,"X","Y"), ordered=TRUE)) %>%
  ggplot(aes(x=pos, y=subclonal_gains)) + 
  geom_area(aes(alpha="sub-clonal", fill="gains")) +
  geom_area(aes(y=clonal_gains, alpha="clonal", fill="gains")) +
  geom_area(aes(y=-subclonal_losses, alpha="sub-clonal", fill="losses")) + 
  geom_area(aes(y=-clonal_losses, alpha="clonal", fill="losses")) + 
  geom_hline(yintercept=0) + 
  geom_vline(data=data.frame(chr=names(centromer_position), pos=centromer_position-chr_centers), 
             aes(xintercept=pos), linetype=3, size=0.5, alpha=0.2) +
  facet_grid(.~chr, scales="free_x", space="free_x") + 
  theme(strip.text.x=element_text(vjust=-102), strip.background.x=element_blank()) + 
  scale_alpha_manual(values=c("clonal"=1, "sub-clonal"=0.45)) +
  scale_fill_manual(values=c("losses"="#3c73a8","gains"="#ec2d01")) + 
  scale_x_continuous(breaks=c(0), labels="") + 
  scale_y_continuous(breaks=seq(-ylm, ylm, by=10), 
                     labels=abs(seq(-ylm, ylm, by=10)), 
                     limits=c(-ylm,ylm)) + 
  xlab("") + ylab("Number of cases") +
  guides(fill="none", alpha="none")

# switch of clipping for the strip texts:
plot_cna_recurence = ggplotGrob(plot_cna_recurence)
for (i in which(grepl("strip-t", plot_cna_recurence$layout$name))){
  plot_cna_recurence$grobs[[i]]$layout$clip = "off"
}

# Plot it!
plot(plot_cna_recurence)

# Save the output
ggsave("results/4_clonality_frequency/genome_wide_clonality_frequency.png", plot_cna_recurence, width=11, height=4.2)

# Output bins to blacklist
trouble_bins = data_per_window_summary[(data_per_window_summary$subclonal_losses > 30 & data_per_window_summary$chr==5) |
                                         (data_per_window_summary$subclonal_losses > 40 & data_per_window_summary$chr==12) |
                                         (data_per_window_summary$subclonal_losses > 40 & data_per_window_summary$chr==19),]
saveRDS(trouble_bins, file = "results/4_clonality_frequency/Bins_to_blacklist_from_recurrence_analysis.rds")
