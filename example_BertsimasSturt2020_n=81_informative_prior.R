  

  ############################################################################
  ### EXAMPLE IN BERTSIMAS AND STURT (2020), N = 81 - WITH INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.05.2021
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  # Revise example in Bertsimas and Sturt (2020)
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )

  
  
  # --------------------------------------------------------------------------
  # Specifications
  # --------------------------------------------------------------------------
  
  n <- 81
  B <- 10^3
  n_boot <- 50
  p_th <- 0.025
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, n - length(z)))
  p_prior <- 1:n / sum(1:n)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  # Ordinary (weighted) bootstrap
  boot_mat_w <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    Boot_w <- boot::boot( data = z,
                          statistic = meanBoot,
                          weights = p_prior,
                          R = B )
    boot_mat_w[ ,j] <- Boot_w$t
  }
  
 
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  bb_mat_w <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = p_prior * n )
    bb_mat_w[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  (toc_bb <- Sys.time() - tic)
  
  
  
  sqrt( mean( apply( boot_mat_w, 2, var ) ) )
  sqrt( mean( apply( bb_mat_w, 2, var ) ) )
  
  
  2
  
  
  
  
  
  
  
  
  
  
  
  
  
  