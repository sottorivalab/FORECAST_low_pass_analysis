# This function runs an ASCAT style algorithm on log2ratio data only
runASCATlp = function(lrrs, fix_ploidy = 2, interval = 0.01, 
                      min_purity = 0.2, max_purity = 1, 
                      max_lrr = Inf, no_fit_psit = 2,
                      pad_ploidy = 0,
                      preset = F,
                      preset_purity,
                      preset_ploidy) {
  
  # This function calculates the distance measure inspired by the ASCAT LogR equation
  # See https://doi.org/10.1073/pnas.1009843107 for the equation
  #
  # rho   = purity of cancer sample
  # psit  = ploidy of cancer cells
  # gamma = 1 for sequencing data
  #
  # Here we calculate sum of squared differences from integer values
  # We also penalise non-positive values (<1)
  fit_lrr = function(lrrs, rho, psit, gamma = 1) {
    
    # Calculate average ploidy of all cells
    psi = (2*(1 - rho)) + (rho*psit)
    
    # Calculate a continuous CN value
    n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
    
    # Calculate what we will compare to
    int_n = abs(round(n))
    
    # Penalise zeros by setting the bottom to 1
    int_n[which(int_n==0)] = 1
    
    # Calculate squared difference
    fit = (n - int_n) ^ 2
    
    # Sum it, that is our 'distance'
    fit = sum(fit) / length(lrrs)
    
    return(fit)
    
  }
  gamma  = 1
  
  # Do we want to remove very high lrrs?
  lrrs_for_fit = lrrs[lrrs < max_lrr]
  
  if(!preset) {
  
  # The range of purities and ploidies we want to search
  rhos   = seq(min_purity, max_purity+interval, by = interval)
  psits  = seq(fix_ploidy - pad_ploidy, fix_ploidy + pad_ploidy, by = interval)
  
  # Run across them all and test their fit
  fits = lapply(rhos, function(r) {
    
    psit_fits =  lapply(psits, function(p) {
      
      fit = fit_lrr(lrrs_for_fit, rho = r, psit = p)
      
    })
    
  })
  
  # Create a matrix from the result
  fit_mat = do.call(rbind, lapply(fits, function(r) {unlist(r)}))
  
  # Name the columns of the fit matrix, to visual fits, make a heatmap of this
  colnames(fit_mat) = psits
  rownames(fit_mat) = rhos
  
  # A function for finding local minima - a 3x3 matrix in which the centre is the minima
  find_local_minima = function(mat) {
    
    # Create the object for storing the results
    res = NULL
    
    # The column numbers to run across, we don't run across the first and last
    for(c in 2:(ncol(mat)-1)) {
      
      # The row numbers to run across, we don't run across the first and last
      for(r in 2:(nrow(mat)-1)) {
        
        # Create a temporary matrix to test for local minima
        test_mat = mat[(r-1):(r+1),(c-1):(c+1)]
        
        # Which entry in the matrix is equal to the minimum in the matrix
        m = which(test_mat==min(test_mat))
        
        # If this is the centre of the matrix and of length one, we have a local minima
        if(m[1]==5 & length(m)==1) {
          
          # Record the psi, purity and value
          hit = c(colnames(test_mat)[2], rownames(test_mat)[2], test_mat[2,2])
          
          # Add to the results object
          res = rbind(res, hit)
          
        }
        
      }
      
    }
    
    # Process it as a data frame
    res = data.frame(psit = as.numeric(res[,1]), purity = as.numeric(res[,2]), 
                     dist = round(as.numeric(res[,3]), digits = 5), 
                     stringsAsFactors = FALSE)
    
    # Order by the distance measured
    res = res[order(res$dist),]
    
    return(res)
    
  }
  
  # Alternate version for a single ploidy state
  find_local_minima_single_ploidy = function(mat) {
    
    # Collect results
    res = NULL
    
    # The row numbers to run across, we don't run across the first and last
    for(r in 2:(nrow(mat)-1)) {
      
      # Create a temporary matrix to test for local minima
      test_mat = mat[(r-1):(r+1),1]
      
      # Which entry in the matrix is equal to the minimum in the matrix
      m = which(test_mat==min(test_mat))
      
      # If this is the centre of the matrix and of length one, we have a local minima
      if(m[1]==2 & length(m)==1) {
        
        # Record the psi, purity and value
        hit = c(psits, names(test_mat)[2], test_mat[2])
        
        # Add to the results object
        res = rbind(res, hit)
        
      }
      
    }
    
    # Process it as a data frame
    res = data.frame(psit = as.numeric(res[,1]), purity = as.numeric(res[,2]), 
                     dist = round(as.numeric(res[,3]), digits = 5), 
                     stringsAsFactors = FALSE)
    
    # Order by the distance measured
    res = res[order(res$dist),]
    
    return(res)
    
  }
  
  # If we are doing a search on a single ploidy, just find the minimum distance
  if(ncol(fit_mat)==1) {
    
    # Which cellularity is correct?
    lms = find_local_minima_single_ploidy(fit_mat)
    
  } else {
    
    # Get local minimas
    lms = find_local_minima(fit_mat)
    
  }
  
  # Get best fit
  best_fit = lms[1,]
  
  # Catch samples without a local minima solution
  if(nrow(lms)==0) {
    
    best_fit = data.frame(psit = no_fit_psit, purity = 1, dist = Inf)
    
  }
  
  }
  
  if(preset) {
    best_fit = data.frame(psit = preset_ploidy, purity = preset_purity, dist = NA)
    fit_mat  = NULL
    lms      = NULL
  }
  
  # Get continuous copy number values for best fit
  rho = best_fit$purity
  psi = (2*(1 - rho)) + (rho*best_fit$psit)
  n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
  
  # Make them integers
  n_int = round(n)
  
  # If we get minus states (i.e. small deletions in impure tumours)
  n_int[n_int<0] = 0
  
  # Return items as a list
  output = list()
  
  # Collect
  output$CN     = n_int
  output$Purity = rho
  output$Psi    = psi
  output$PsiT   = best_fit$psit
  output$contCN = n
  output$fit    = fit_mat
  output$allsol = lms
  
  return(output)
  
}
