  
  
  ############################################################################
  ### EXAMPLE Cales Chalkis Emiris 2021
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.08.2021
  # First version:    13.08.2021
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
  
  
  
 
  
  
  
  
  
  # Question: Is algo 3 in CalesChalkisEmiris2021 generalizable to weighted moments?
  # Partial answer: We know that the closed form expressions for the first four
  # moments hold when replacing functions mean, var, skew and kurt with their 
  # weighted counterparts (see below).
  
  
  n <- 10
  df <- 4
  n_sim <- 10^3
  n_boot <- 50
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  # --------------------------------------------------------------------------
  # Unweighted
  # --------------------------------------------------------------------------

  # Compute M_2 by Theorem 15
  M2 <- sum( (z - mean(z))^2 ) / (n * (n + 1) )
    
  # k-th moment:
  kvec <- 3:10
  lMk <- lapply( kvec, function(k) { getMk( z = z, k = k ) } )
  mk <- setNames( unlist(lMk), paste0("M", kvec) )
  mk
  lB <- lapply( kvec, function(k) { attr(getMk( z = z, k = k ), "B" ) } )
  B <- setNames( unlist(lB), paste0("M", kvec) )
  
  
  # Compare with closed form solution
  M4FUN( x = z, M2 = M2 )
  M3FUN( x = z, M2 = M2 )
  
  s_unb <- var(z)
  s_unb * ( (2 / (n + 1)) - 1 / n )
  M2
  
  A <- ( (2 / (n + 1)) - 1 / n )
  cbind( mk, 
         1 / M2^(kvec/2) * B, 
         1 / s_unb^(kvec/2) * B / A^(kvec/2),
         sqrt(s_unb)^(-kvec) * B / A^(kvec/2) )  # same same
  
  
  
  
  # --------------------------------------------------------------------------
  # Weighted
  # --------------------------------------------------------------------------
  
  # Parameter for asymmetric Dirichlet
  
  # wghts <- 1:n / sum(1:n)
  wghts <- rep(1/n, n)
  alpha <- wghts * n 

  # Weighted M2
  mu_wt <- sum( wghts * z )
  M2_wt <- sum( wghts * (z - mu_wt)^2 ) / (n + 1)
  
  # Weighted M3
  M3_wt <- (1 / M2_wt^(3/2)) * (2 / ((n + 1) * (n + 2))) * sum( wghts * (z - mu_wt)^3 )
  
  # Weighted M4
  M4_wt <- (1 / M2_wt^2) * 
             (1 / ((n + 1) * (n + 2) * (n + 3))) * 
               ( 6 * sum( wghts * (z - mu_wt)^4 ) + 3 * sum(wghts^2) * sum( (z - mu_wt)^2 )^2 )
  
  M2_wt
  M3_wt
  M4_wt

  
  
  # Weighted Bayesian bootstrap
  bb_mat <- matrix( NA, nrow = n_sim, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = n_sim, alpha = alpha )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  var_bb <- apply( bb_mat, 2, var )
  range( var_bb )
  M2_wt; mean( var_bb )
  
  sk_bb <- apply( bb_mat, 2, skewFUN )
  range( sk_bb )
  M3_wt; mean( sk_bb )
  
  kurt_bb <- apply( bb_mat, 2, kurtFUN )
  range( kurt_bb )
  M4_wt; mean( kurt_bb )
  
  # m5_bb <- apply( bb_mat, 2, fifthMomFUN )
  # range( m5_bb )
  # M5_wt; mean( m5_bb )
  
  
  
  Mk3 <- getMk( z = z, k = 3 )
  Mk3
  M3FUN( x = z )  # same same
  M3_wt
  
  attr(Mk3, "S")
  attr(Mk3, "B")
  
  1 / M2^(3/2) * attr(Mk3, "B"); M3FUN( x = z )

  num <- choose(n = n - 1 + 3, 3)
  1 / M2_wt^(3/2) * attr(Mk3, "S") / num; M3_wt
  
  M3_wt / (1 / M2_wt^(3/2)) * num
  attr(Mk3, "S")
  
  
  
  
  Mk4 <- getMk( z = z, k = 4 )
  Mk4
  M4FUN( x = z )  # same same
  mean(kurt_bb)
  M4_wt
  
  
  Mk5 <- getMk( z = z, k = 5 )
  Mk5
  m5_bb <- apply( bb_mat, 2, momFUN, k = 5 )
  mean( m5_bb )
  range( m5_bb )
  
  
  k <- 3; as.numeric( getMk( z = z, k = k ) ); mean (apply( bb_mat, 2, momFUN, k = k ) )
  k <- 4; as.numeric( getMk( z = z, k = k ) ); mean (apply( bb_mat, 2, momFUN, k = k ) )
  k <- 5; as.numeric( getMk( z = z, k = k ) ); mean (apply( bb_mat, 2, momFUN, k = k ) )
  k <- 6; as.numeric( getMk( z = z, k = k ) ); mean (apply( bb_mat, 2, momFUN, k = k ) )
  k <- 7; as.numeric( getMk( z = z, k = k ) ); mean (apply( bb_mat, 2, momFUN, k = k ) )
 
  
  
  
  attr(Mk3, "M2")
  S <- attr(Mk3, "S")
  S / attr(Mk3, "B")
  choose(n = n - 1 + k, k)
  choose(n = n, k)
  
  
  1 / attr(Mk3, "M2")^(3/2) * S / choose(n = n - 1 + k, k)
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # xxx
  # --------------------------------------------------------------------------
  
  n <- 8
  B <- 10^3
  n_boot <- 50
  df <- 5
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  
  s_unb <- var(z)
  s_ml <- sum( (z - mean(z))^2 ) / n
  sk <- skewFUN( x = z )
  ks <- kurtFUN( x = z )
  
  # Bayesian bootstrap
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = rep(1, n) )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  
  sk
  M3FUN( x = z, M2 = NULL )
  
  ks
  M4FUN( x = z, M2 = NULL )
  14 / 3
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BB with non-informative and informative prior, exact first four moments
  # Loop over n
  # --------------------------------------------------------------------------
  

  n_sim <- 10^4 + 1
  n_boot <- 50
  n_vec <- seq(from = 2, to = 83, by = 1)
  
  for ( n in n_vec ) {
    
    
    # --------------------------------------------------------------------------
    # Data
    # --------------------------------------------------------------------------
    
    df <- 4
    set.seed(1111)
    x <- rt( n = n, df = df )
    z <- x^2
    
    
    # --------------------------------------------------------------------------
    # Bayesian bootstrap
    # --------------------------------------------------------------------------
    
    # Non-informative
    
    bb_mat <- matrix( NA, nrow = n_sim, ncol = n_boot )
    lP <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = n_sim, alpha = rep(1, n) )
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    # Informative prior

    # Parameter for asymmetric Dirichlet
    wghts <- 1:n / sum(1:n)
    alpha <- wghts * n 
    
    bb_mat_info <- matrix( NA, nrow = n_sim, ncol = n_boot )
    lP <- list()
    lNorm <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = n_sim, alpha = alpha )
      lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
      bb_mat_info[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    # Expected 2-norm
    expected_2normSq_bb <- unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) )
    
    
    
    # --------------------------------------------------------------------------
    # Variance
    # --------------------------------------------------------------------------
    
    s_unb <- var(z)
    s_unb_w <- as.numeric( cov.wt( x = as.matrix(z, ncol = 1), 
                                  wt = wghts,
                                  cor = FALSE,
                                  center = TRUE,
                                  method = "unbiased" )$cov )
    s_ml <- sum( (z - mean(z))^2 ) / n
    s_ml_w <- as.numeric( cov.wt( x = as.matrix(z, ncol = 1), 
                                  wt = wghts,
                                  cor = FALSE,
                                  center = TRUE,
                                  method = "ML" )$cov )
                        
    # Bayesian bootstrap (BB)
    stderr_bb <- apply( bb_mat, 2, sd )
    
    # Exact BB
    stderr_exact_bb <- sqrt( s_unb * ( (2 / (n + 1)) - 1 / n ) )
    # stderr_exact_bb <- sqrt( s_ml * ( 1 + (2 / (n + 1) - 1) / (1 - 1/n) ) )
    stderr_exact_cales <- sqrt( sum((z - mean(z))^2) / (n * (n + 1)) )  # same same
    stderr_exact_cales <- sqrt( var(z) * ((n - 1) / n) / (n + 1) )  # same same
    
    sqrt( s_ml / (n + 1) )
    stderr_exact_bb
    
    
    
    # Bayesian bootstrap (BB) - Informative prior
    stderr_bb_info <- apply( bb_mat_info, 2, sd)
    
    # Exact BB
    betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
    beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
    # beta_moments
    # (alpha^2 + alpha) / (n^2 + n) # same same
    # stderr_exact_bb_info <- sqrt( s_unb_w * ( sum(beta_moments) - lpNorm(wghts, p = 2)^2 ) )
    stderr_exact_bb_info <- sqrt( s_ml_w * ( 1 + (sum(beta_moments) - 1) / (1 - lpNorm(wghts, p = 2)^2) ) )
    # stderr_exact_bb_info <- sqrt( s_ml_w / (n + 1) )

    stderr_exact_bb_info
    sqrt( s_ml_w / (n + 1) )  # same same
    mean( stderr_bb_info )  # close enough
    

    
    
    
    # --------------------------------------------------------------------------
    # Skewness
    # --------------------------------------------------------------------------
    
    # Bayesian bootstrap (BB)
    M3_bb <- apply( bb_mat, 2, FUN = skewFUN )
    
    # Exact BB
    M3_exact_bb <- M3FUN( x = z, M2 = stderr_exact_bb^2 )
    M3_exact_bb2 <- M3FUN( x = z, M2 = NULL )
    M3_exact_bb; M3_exact_bb2 # same same
    mean(M3_bb)
    
    
    # Informative prior
    M3_bb_info <- apply( bb_mat_info, 2, FUN = skewFUN )
    
    
    sk_w <- skew.wt( x = z, wt = wghts )
    sk_w * (2 * sqrt(n + 1)) / (n + 2)
    mean(M3_bb_info)  # close enough
    M2_w <- stderr_exact_bb_info^2
    mu_w <- sum(wghts * z)
    (1 / M2_w^(3/2)) * (2 / ((n + 1) * (n + 2))) * sum( wghts * (z - mu_w)^3 ) # same same
    
    
    
    # debugonce( getMk )
    Mk <- getMk( z = z, k = 3, M2 = M2_w )
    B <- attr(Mk, "B")
    1 / M2_w^(3/2) * B
    Mk
    
    
    
    
    # --------------------------------------------------------------------------
    # Kurtosis
    # --------------------------------------------------------------------------
    
    # Bayesian bootstrap (BB)
    M4_bb <- apply( bb_mat, 2, FUN = kurtFUN )
    
    # Exact BB
    M4_exact_bb <- M4FUN( x = z, M2 = stderr_exact_bb^2 )
    M4_exact_bb2 <- M4FUN( x = z, M2 = NULL )
    M4_exact_bb; M4_exact_bb2 # same same
    mean(M4_bb)
    
    
    # Informative prior
    M4_bb_info <- apply( bb_mat_info, 2, FUN = kurtFUN )
    
    
    ks_w <- kurt.wt( x = z, wt = wghts )
    (2 * ks_w + n) * (3 * (n + 1)) / ((n + 2) * (n + 3))
    mean(M4_bb_info) # close enough
    M2_w <- stderr_exact_bb_info^2
    mu_w <- sum(wghts * z)
    (1 / M2_wt^2) * 
      (1 / ((n + 1) * (n + 2) * (n + 3))) * 
      ( 6 * sum( wghts * (z - mu_wt)^4 ) + 3 * sum(wghts^2) * sum( (z - mu_wt)^2 )^2 )

    
    
    
    
    
  }
  
  
  
  
  
  
    
    
    
    
    
    
    