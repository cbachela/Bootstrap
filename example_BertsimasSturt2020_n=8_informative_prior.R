  

  ############################################################################
  ### EXAMPLE IN BERTSIMAS AND STURT (2020), N = 8, INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.07.2021
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  # Show that the example in Bertsimas and Sturt (2020) can be formulated as 
  # a n = 8 dimensional problem with non-uniform probabilities.
  
  # NO, it can not!
  # Standard deviation is much higher in the n=8 dimensional example and quantile
  # cannot be computed exactly because the exact distribution is a step-function.
  
  
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
  
  n <- 8
  B <- 10^3
  n_boot <- 10^3
  p_th <- 0.025
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, 1)
  p_prior <- c(rep(1/81, 7), 74/81)
 
  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
 
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )

  # Reorder results by decreasing total probability
  ordering <- rev( order( prob_tot) )
  parts_mat <- parts_mat[ ,ordering]
  draws <- draws[ ,ordering]
  prob <- prob[ordering]
  prob_tot <- prob_tot[ordering]
  prob_dmult <- prob_dmult[ordering]

  # Euclidean norm of partitions
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )

  
  
  
  
  # Enumerate all permutations within each level set
  
  FUN <- function(i)
  {
    P <- RP::permutations(parts_mat[ ,i])
    P <- matrix( P[!duplicated(P), ], ncol = n )
  }
  tic <- Sys.time()
  lP <- lapply( 1:ncol(parts_mat), FUN = FUN )
  (toc <- Sys.time() - tic)
  
  # Probability of each grid node, grouped by partitions
  lprob <- lapply( lP, function(P) { apply( P, 1, dmultinom, prob = p_prior, log = FALSE ) } )
  
  # Theta
  lTheta <- lapply( lP, FUN = function(P) { (P / n) %*% z } )
  
  # Variance (biased) of theta's by partition
  theta_var <- unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(theta))^2 ) / length(theta) } ) )
  
 
  
  
  # Exact expected 2-norm under prior probabilities
  lp_lP <- unlist( lapply( lP, function(P) { apply( P / n, 1, lpNorm, p = 2 ) } ) )
  elp2 <- sum( lp_lP^2 * unlist(lprob) )
  elp2
  
  
  
  # Ordinary (weighted) bootstrap
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_w <- boot_mat
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    Boot_w <- boot::boot( data = z,
                        statistic = meanBoot,
                        weights = p_prior,
                        R = B )
    boot_mat_w[ ,j] <- Boot_w$t
    
  }
  
  # Unweighted
  mean( apply( boot_mat, 2, var ) )
  ( sum( lp_parts^2 * prob_tot ) - 1/n ) * var(z) # close enough
  
  # Weighted
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                    wt = p_prior, 
                    cor = FALSE, 
                    center = TRUE, 
                    method = "unbiased" )
  Vz_unb$cov * ( elp2 - lpNorm(p_prior, p = 2)^2 )
  mean( apply( boot_mat_w, 2, var ) ) # close enough
  
  

  
  
  # # Exact quantile -  find the quantile H(0.025) iteratively
  # theta_vec <- unlist( lTheta )
  # p_vec <- unlist( lprob )
  # err <- 1
  # tol <- 1e-12
  # p_th <- 0.025
  # z0 <- quantile(z, p_th)
  # z_max <- max(z)
  # i <- 1
  # tic <- Sys.time()
  # while( err > tol ) {
  #   idx <- which( theta_vec <= z0 )
  #   p <- sum(p_vec[idx])
  #   err <- abs(p - p_th)
  #   if ( p > p_th ) {
  #     z0 <- z0 - 1/(2^i) * z_max
  #   } else {
  #     z0 <- z0 + 1/(2^i) * z_max
  #   }
  #   i <- i + 1
  # }
  # (toc <- Sys.time() - tic)
  # i
  # 
  # 
  # 
  # # Compare quantile values
  # z0; p
  # quantile( Boot_w$t, p_th )  
  # mean( apply(boot_mat_w, 2, quantile, p_th) )
  # # True value is 37.08
  
  
  
  
  
  
  # Compare to original n=81 example
  env <-  readRDS( file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  ls(env)  
    
  
  
  
  sqrt( Vz_unb$cov * ( elp2 - lpNorm(p_prior, p = 2)^2 ) )
  sqrt( mean( apply( boot_mat_w, 2, var ) ) )
  sqrt( sum( env$prob_tot * (env$lp_parts^2 * var(env$z) - min(env$lp_parts^2) * var(env$z)) ) )
  
  # different!
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian approach
  # --------------------------------------------------------------------------
  
  # We compare two procedures:
  # 1) BB based on n = 81
  # 2) BB based on n = 8 with informative prior

  
  n_boot <- 50

  # 1) BB based on n = 81
  
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = p_prior * n )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  sqrt( mean( apply( bb_mat, 2, var ) ) ) # different!
  env$SE   
  
  
  
  
  
  
  
  
  
  
  
  
  