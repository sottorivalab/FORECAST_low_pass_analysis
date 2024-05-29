library(phytools)

# Get patients
medicc_patients = gsub("_regions.tsv", "", list.files("~/remotes/medicc_trees/forecast", pattern = "_regions.tsv"))

# Run through them
phylosig = lapply(medicc_patients, function(patient) {
  
  # Read in the file
  tree = read.tree(paste0("~/remotes/medicc_trees/forecast/",patient,"_regions_final_tree.new"))
  
  # Make a short df of tree samples
  sample_side = data.frame(Sample = tree$tip.label, stringsAsFactors = F)
  
  # Get sidedness information
  sample_side$Region = unlist(lapply(strsplit(sample_side$Sample, split = "_"), function(i) unlist(strsplit(i[2], split = ""))[1]))
  sample_side$Region_score = ifelse(sample_side$Region=="R", 1, ifelse(sample_side$Region=="L", yes = 0, no = 0.5))
  
  print(patient)
  if(!any(sample_side$Region[!is.na(sample_side$Region)] %in% c("R", "L"))) {print("non L R sample")}
  
  pheno = sample_side$Region_score
  names(pheno) = sample_side$Sample
  
  #tree$edge.length[length(tree$edge.length)] = 1
  
  res = try(phylosig(tree, pheno, method='lambda', test=T))
  
  return(res)
  
})

results = data.frame(Patient = medicc_patients, 
                     p = signif(unlist(lapply(phylosig, function(i) {ifelse(class(i)=="try-error", NA, i$P)}))), 
                     lambda = signif(unlist(lapply(phylosig, function(i) {ifelse(class(i)=="try-error", NA, i$lambda)}))))

# Output for manuscript and survival analysis
options(scipen = 999)
write.csv(results, file = "results/8_collect_metrics/FORECAST_phylosig_sidedness_analysis.csv", quote = F, row.names = F)
