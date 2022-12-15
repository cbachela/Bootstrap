  
  
  ############################################################################
  ### NUMERICAL EXAMPLE OF BERTSIMAS AND STURT 2020
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.09.2022
  # First version:    16.06.2022
  # --------------------------------------------------------------------------
  
  # Replicate example in Bertsimas and Sturt (2020)
  
  
 
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(partitions)
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
  n_boot <- 10^2
  p_th <- 0.025
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, n - length(z)))
  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
  # parts_mat <- partitions::parts( n )
  # draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  # prob <- draws2prob( draws = draws )
  # prob_tot <- prob * draws[2, ]
  # prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  # 
  # # Reorder results by decreasing total probability
  # ordering <- rev( order( prob_tot) ) 
  # parts_mat <- parts_mat[ ,ordering]
  # draws <- draws[ ,ordering]
  # prob <- prob[ordering]
  # prob_tot <- prob_tot[ordering]
  # prob_dmult <- prob_dmult[ordering]
  # 
  # # Euclidean norm of partitions
  # lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  # 
  # 
  # head( cbind(t(draws), log(prob), log(prob_dmult), prob_tot) )
  #
  #
  # env <- new.env()
  # 
  # env$n <- n
  # env$z <- z
  # env$B <- B
  # env$n_boot <- n_boot
  # env$p_th <- p_th
  # env$parts_mat <- parts_mat
  # env$draws <- draws
  # env$prob <- prob
  # env$prob_tot <- prob_tot
  # enb$lp_parts <- lp_parts
  # 
  # saveRDS( env, file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  
  
  
  env <- readRDS( file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  n <- env$n
  z <- env$z
  B <- env$B
  p_th <- env$p_th
  n_boot <- env$n_boot
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  prob_tot <- env$prob_tot
  lp_parts <- env$lp_parts
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = rep(1, n) )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  (toc_bb <- Sys.time() - tic)
  
  
  
  
  
  
  # Save
  env$boot_mat <- boot_mat
  env$bb_mat <- bb_mat
  saveRDS( env, file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  
  
  
  
  
  