  
  
  ############################################################################
  ### EXACT STANDARD ERROR - INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.06.2021
  # First version:    16.06.2021
  # --------------------------------------------------------------------------
  
  
  # Show that the standard error of the Bayesian bootstrap with informative
  # prior can be computed exactly.
  
  

  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  df <- 10
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  m <- n

  # Weights
  wghts <- 1:n / sum(1:n)
  
  # Dirichlet parameter
  alpha <- wghts * n 
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary weighted bootstrap
  # --------------------------------------------------------------------------
  
 
  tic <- Sys.time()
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B,
                        weights = wghts )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  
  # Alternatively:
  tic <- Sys.time()
  boot_mat2 <- matrix( NA, nrow = B, ncol = n_boot )
  lNorm_boot <- list()
  for ( j in 1:n_boot ) {
    gridnodes <- matrix( 0, nrow = B, ncol = n )
    for ( i in 1:B ) {
      idx <- sample( 1:n, size = m, replace = TRUE, prob = wghts )
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / m
    }
    lNorm_boot[[j]] <- apply( gridnodes, 1, lpNorm, p = 2 )
    boot_mat2[ ,j] <- apply( gridnodes, 1, function(p) { sum(p * z) } )
  }
  (toc_boot2 <- Sys.time() - tic)
  

  
  # Expected 2-norm
  expected_2normSq_boot <- unlist( lapply( lNorm_boot, FUN = function(x) { mean(x^2) } ) )
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian asymmetric bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  (toc_bb <- Sys.time() - tic)
  
  
  # Expected 2-norm
  expected_2normSq_bb <- unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) )
  
  
  
  # --------------------------------------------------------------------------
  # Standard error
  # --------------------------------------------------------------------------
  
  # Weighted sample variance (unbiased)
  Vz <- cov.wt( x = matrix(z, ncol = 1), 
                wt = wghts, 
                cor = FALSE, 
                center = TRUE, 
                method = "unbiased" )
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
    
  # Ordinary bootstrap
  stderr_boot <- apply( boot_mat, 2, sd )
  stderr_boot2 <- apply( boot_mat2, 2, sd )
  
  # Exact ordinary bootstrap 
  # Actually, pseudo exact bcs expected norm is only estimated!
  stderr_parts <- sqrt( Vz$cov * ( mean(expected_2normSq_boot) - lp_wbar^2 ) )
  stderr_parts2 <- sqrt( Vz$cov * ( mean(expected_2normSq_boot) - lp_wbar^2 ) )
  
  Vz_ml$cov * (1 + (mean(expected_2normSq_boot) - 1) / (1 - lp_wbar^2))
  stderr_parts^2
  
 
  # Bayesian bootstrap (BB)
  stderr_bb <- apply( bb_mat, 2, sd )
 
  # Exact BB
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  # beta_moments
  # (alpha^2 + alpha) / (n^2 + n) # same same
  stderr_exact_bb <- sqrt( Vz$cov * ( sum(beta_moments) - lpNorm(wghts, p = 2)^2 ) )

  
  
  
  # Compare
  SE <- setNames( c(mean(stderr_boot), 
                    mean(stderr_boot2),
                    stderr_parts, 
                    stderr_parts2,
                    mean(stderr_bb), 
                    stderr_exact_bb,
                    sqrt(Vz$cov) / sqrt(n)),
                  c("avg_boot",
                    "avg_boot2",
                    "exact_boot", 
                    "exact_boot2",
                    "avg_bb", 
                    "exact_bb",
                    "normal") )
  SE
  
  
  
  
  # # --------------------------------------------------------------------------
  # # Skewness
  # # --------------------------------------------------------------------------
  # 
  # # Ordinary bootstrap
  # M3_boot <- apply( boot_mat, 2, skewFUN )
  # M3_boot2 <- apply( boot_mat2, 2, skewFUN )
  # 
  # # Exact ordinary bootstrap
  # M3_parts <- M3FUN( x = z, M2 = stderr_parts^2 )
  # M3_parts2 <- M3FUN( x = z, M2 = stderr_parts2^2 )
  # 
  # # Bayesian bootstrap (BB)
  # M3_bb <- apply( bb_mat, 2, FUN = skewFUN )
  # 
  # # Exact BB
  # M3_exact_bb <- M3FUN( x = z, M2 = stderr_exact_bb^2 )
  # 
  # 
  # # Compare
  # M3 <- setNames( c(mean(M3_boot), 
  #                   mean(M3_boot2),
  #                   M3_parts, 
  #                   M3_parts2,
  #                   mean(M3_bb), 
  #                   M3_exact_bb,
  #                   0),
  #                 c("avg_boot",
  #                   "avg_boot2",
  #                   "exact_boot", 
  #                   "exact_boot2",
  #                   "avg_bb", 
  #                   "exact_bb",
  #                   "normal") )
  # M3
  # 
  # 
  # 
  # 
  # 
  # # --------------------------------------------------------------------------
  # # Kurtosis
  # # --------------------------------------------------------------------------
  # 
  # # Ordinary bootstrap
  # M4_boot <- apply( boot_mat, 2, kurtFUN )
  # M4_boot2 <- apply( boot_mat2, 2, kurtFUN )
  # 
  # # Exact ordinary bootstrap
  # M4_parts <- M4FUN( x = z, M2 = stderr_parts^2 )
  # M4_parts2 <- M4FUN( x = z, M2 = stderr_parts2^2 )
  # 
  # # Bayesian bootstrap (BB)
  # M4_bb <- apply( bb_mat, 2, FUN = kurtFUN )
  # 
  # # Exact BB
  # M4_exact_bb <- M4FUN( x = z, M2 = stderr_exact_bb^2 )
  # 
  # 
  # # Compare
  # M4 <- setNames( c(mean(M4_boot), 
  #                   mean(M4_boot2),
  #                   M4_parts, 
  #                   M4_parts2,
  #                   mean(M4_bb), 
  #                   M4_exact_bb,
  #                   # kurtFUN( rnorm(10^7, mean(z), sd(z) / sqrt(n)))),
  #                   3),
  #                 c("avg_boot",
  #                   "avg_boot2",
  #                   "exact_boot", 
  #                   "exact_boot2",
  #                   "avg_bb", 
  #                   "exact_bb",
  #                   "normal") )
  # M4
  # 
  # 
  # SE
  # M3
  # M4
  # 
  # (3 * df - 6) / (df - 4) # Analytic kurtosis of t-dist
  
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  env <- new.env()
  
  env$n <- n
  # env$m <- m
  env$B <- B
  env$n_boot <- n_boot
  # env$lp_parts <- lp_parts
  # env$prob_tot <- prob_tot
  env$avg_expected_2normSq_boot <- mean(expected_2normSq_boot)
  env$avg_expected_2normSq_bb <- mean(expected_2normSq_bb)
  env$sum_beta_moments <- sum(beta_moments)
  env$SE <- SE
  # env$M3 <- M3
  # env$M4 <- M4
  
  
  saveRDS( object = env, file = paste0(wd, "waRehouse/exact_std_err_info_prior_n=", n) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    n_vec <- seq(from = 2, to = 415, by = 1)
    lSE <- list()
    e2n <- numeric(length(n_vec))  # expected 2-norm
    a2nboot <- numeric(length(n_vec))  # average 2-norm ordinary bootstrap
    a2nbb <- numeric(length(n_vec))  # average 2-norm bayesian bootstrap
    for ( i in seq(along = n_vec) ) {
      env <- readRDS( file = paste0(wd, "waRehouse/exact_std_err_info_prior_n=", n_vec[i]) )
      lSE[[i]] <- env$SE
      e2n[i] <- env$sum_beta_moments
      a2nboot[i] <- env$avg_expected_2normSq_boot
      a2nbb[i] <- env$avg_expected_2normSq_bb
    }
   
    
    lp_wbar <- unlist( lapply( n_vec, FUN = function(n) { lpNorm( 1:n / sum(1:n), p = 2 ) } ) )
    e2n_boot <- unlist( lapply( n_vec, FUN = function(n) { 1/n - lp_wbar[n-1]^2/n + lp_wbar[n-1]^2 } ) )
    A <- cbind( uninf = 1 + (2/n_vec - 1/n_vec^2 - 1) / (1 - 1/n_vec),
                inf_approx = 1 + (a2nboot - 1) / (1 - lp_wbar^2),
                inf_exact = 1 + (e2n_boot - 1) / (1 - lp_wbar^2),
                xxx = 1 / n_vec )
    head(A)
    plot( as.timeSeries(A), plot.type = "single" )
    
    cbind( a2nboot, e2n_boot)
    
    
    
    two_norm <- cbind( expected = e2n, 
                       average_boot = a2nboot, 
                       exact_boot = e2n_boot,
                       average_bb = a2nbb )
    colors <- 1:ncol(two_norm)
    plot( as.timeSeries(two_norm), plot.type = "single", col = colors )
    legend( "topright", colnames(two_norm), lwd = 2, col = colors, text.col = colors, bty = "n" )
    
    head(two_norm)
    
    
    
    SE <- t( do.call( cbind, lSE ) )
    rownames(SE) <- paste0("n=", n_vec[1:nrow(SE)])
    
    FUN <- function(n) 
    {
      s_tilde <- cov.wt( x = matrix(z[1:n], ncol = 1), 
                         wt = 1:n / sum(1:n),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML" )$cov
      ans <- sqrt( as.numeric(s_tilde) * A[n-1, ] )
      return( ans )
    }
    SE <- cbind( SE, do.call( rbind, lapply( n_vec, FUN = FUN ) ) ) 
    
    SE
    head(SE)
    tail(SE)
    
    
    
    colors <- c("red", "darkgreen", "black")
    plot( x = n_vec[1:nrow(SE)], y = SE[ ,1], ylim = range(SE), type = "o" )
    for ( i in 1:ncol(SE) ) {
      lines( x = n_vec[1:nrow(SE)], y = SE[ ,i], col = colors[i] ) 
    }
    
    plot( x = n_vec[1:nrow(SE)], y = (SE[ ,"avg_boot"] - SE[ ,"exact_boot"]), type = "o" )
    abline( h = 0, col = "grey" )

        plot( x = n_vec[1:nrow(SE)], y = (SE[ ,"avg_bb"] - SE[ ,"exact_bb"]), type = "o" )
    abline( h = 0, col = "grey" )
    
    plot( x = n_vec[1:nrow(SE)], y = (SE[ ,"normal"] - SE[ ,"exact_bb"]), type = "o" )
    abline( h = 0, col = "grey" )  
    
    
    
    
    
    
  }
  
  
  
  
 
  
  
  
  
  
  
  
  
  