  
  
  ############################################################################
  ### EXACT STANDARD ERROR - RUNFILE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.05.2021
  # First version:    22.05.2021
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  
  
  # Run source script
  B <- 10^4 + 1
  n_boot <- 50
  n_vec <- seq(from = 2, to = 83, by = 1)
  # n_vec <- seq(from = 2, to = 500, by = 1)

  for ( n in n_vec ) {
    
      m <- ceiling( n * 0.8 )
      source( paste0(wd, "exact_standard_error.R") )
      # source( paste0(wd, "exact_standard_error_informative_prior.R") )
    
  }
  
  
  
  
  
  
  
  
  
  
  