  
  
  ############################################################################
  ### EXACT MOMENTS - BAYESIAN BOOTSTRAP - INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.04.2022
  # First version:    20.10.2021
  # --------------------------------------------------------------------------
  

  # Sections:
  #
  # Case 0: Symmetric alpha and sum(alpha) = n
  # Case 1: Asymmetric alpha, but sum(alpha) = n
  # Case 2a: Symmetric alpha, but sum(alpha) > n
  # Case 2b: Symmetric alpha with sum(alpha) < n 
  # Case 3a: Asymmetric alpha, but sum(alpha) > n
  # Case 3b: Asymmetric alpha with sum(alpha) < n 
  # Case 4a: FEV-Bias, symmetric
  # Case 4b: FEV-Bias, Asymmetric
 
  
  
  
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
  
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 15
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1:10)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  # Findings:
  
  # The exact formula 
  # \mu_2(\theta) = m_2 ( \mathds{E}(||\omega||_m^m) - ||\bar{\omega}||_m^m )
  # holds for all cases 0-3.
  
  # Analytic expressions for M2-M4 fail to hold when \sum(alpha) <> n.
  # This is because the closed form expressions implicitely assume 
  # \mathds{E}(||\omega||_m^m) according to \sum_{i=1}^n \alpha_i = n.
  
  
  
  
  
  
  # The closed-form equation for the variance holds in for any parametrization of the Dirichlet.
  
  
  # Case 1, asymmetric alpha with sum(alpha) = 1: all four (scaled and centred) moments
  # of the sample mean are correctly computable with the weighted moment functions (i.e.,
  # scaling factors are unchanged).
  
  # Case 2 the closed-form variance of the mean computations are correct when using the 
  # expected squared 2-norms. Not so with the s(w)/(n+1) rule because the distances term 
  # is not constant.
  # I.e., the scaling factor (1 + (sum(beta_moments) - 1) / (1 - lp_wbar^2) ) != 1 / (n + 1)
  # --> SE can be computed exactly when using the 2-norm based approach.
    
  # This also seems to be true for case 3 and 4a. Not so for the asymmetric FEV case 4b.
  
  
  
  # Question: 
  # Can we find the correct scaling factor for cases 2-4a s.t. we can adjust functions
  # M3BBFUN and M4BBFUN? --> Actually, can we find correct adjustment for function
  # getMk?
  
  
 
  
  # --------------------------------------------------------------------------
  # Case 0: Symmetric alpha and sum(alpha) = n
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  wghts <- rep(1/n, n)
  alpha <- wghts * n 
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  lNorm3 <- lNorm4 <- lNorm5 <- lNorm6 <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    lNorm3[[j]] <- apply( lP[[j]], 1, lpNorm, p = 3 )
    lNorm4[[j]] <- apply( lP[[j]], 1, lpNorm, p = 4 )
    lNorm5[[j]] <- apply( lP[[j]], 1, lpNorm, p = 5 )
    lNorm6[[j]] <- apply( lP[[j]], 1, lpNorm, p = 6 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected 2-norms
  expected_2normSq <- mean( unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) ) )
  expected_3norm3 <- mean( unlist( lapply( lNorm4, FUN = function(x) { mean(x^3) } ) ) )
  expected_4norm4 <- mean( unlist( lapply( lNorm4, FUN = function(x) { mean(x^4) } ) ) )
  expected_5norm5 <- mean( unlist( lapply( lNorm4, FUN = function(x) { mean(x^5) } ) ) )
  expected_6norm6 <- mean( unlist( lapply( lNorm4, FUN = function(x) { mean(x^6) } ) ) )
  
  # Exact expected (2-norm)^2, (3-norm)^3, (4-norm)^4
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
  
  # Variance of sample mean
  M2 <- M2Exact( z = z, w_bar = wghts, exp2norm2 = sum(beta_moments_2) )
  M2
  M2Analytic( z = z, w_bar = wghts, method = "Bayesian" )
  Vz_ml$cov / (n + 1)
  mean( apply( bb_mat, 2, var ) )
  
  
  # Skewness of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^3 ) } ) )
  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  M3Exact( z = z, 
           w_bar = wghts, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2) )
  # M3BBFUN( x = z, M2 = M2, wghts = wghts )
  # M3BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, skewFUN ) )
  
  
  # Kurtosis of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^4 ) } ) )
  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  # M4BBFUN( x = z, M2 = M2, wghts = wghts )
  # M4BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, kurtFUN ) )  
  

  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 1: Asymmetric alpha, but sum(alpha) = n
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  wghts <- 1:n / sum(1:n)
  w_bar <- wghts
  alpha <- wghts * n 
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  lNorm3 <- lNorm4 <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    lNorm3[[j]] <- apply( lP[[j]], 1, lpNorm, p = 3 )
    lNorm4[[j]] <- apply( lP[[j]], 1, lpNorm, p = 4 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Exact expected (2-norm)^2, (3-norm)^3, (4-norm)^4
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  sum(beta_moments_2); mean(unlist(lNorm)^2)
  sum(beta_moments_3); mean(unlist(lNorm3)^3)
  sum(beta_moments_4); mean(unlist(lNorm4)^4)
  
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
 
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
  
  # Variance of sample mean
  M2 <- M2Exact( z = z, w_bar = wghts, exp2norm2 = sum(beta_moments_2) )
  M2
  M2Analytic( z = z, w_bar = wghts, method = "Bayesian" )
  Vz_ml$cov / (n + 1)
  mean( apply( bb_mat, 2, var ) )
  
  
  # Skewness of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^3 ) } ) )
  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  M3Exact( z = z, 
           w_bar = wghts, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2) )
  # M3BBFUN( x = z, M2 = M2, wghts = wghts )
  # M3BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, skewFUN ) )
  
  
  # Kurtosis of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^4 ) } ) )
  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  # M4BBFUN( x = z, M2 = M2, wghts = wghts )
  # M4BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, kurtFUN ) )  
  
  
 
  
  
  
  
  

  # --------------------------------------------------------------------------
  # Case 2: Symmetric alpha, but sum(alpha) <> n 
  #         I.e., unweighted problem but different weight concentration (hence
  #         very similar to case 4a except that E(||w||_p^p) is known and 
  #         w \sim Dirichlet(\alpha)
  # --------------------------------------------------------------------------
 
  # --------------------------------------------------------------------------
  # Case 2a:  Symmetric alpha, but sum(alpha) > n 
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  alpha <- rep(13, n)
  wghts <- alpha / sum(alpha)
  w_bar <- wghts
    
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected squared 2-norm
  expected_2normSq <- mean( unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) ) )
  
  # Exact expected (2-norm)^2, (3-norm)^3, (4-norm)^4
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  Vz_ml$cov; sum( wghts * (z - sum(wghts*z))^2 )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
 
  # Variance of sample mean
  M2 <- M2Exact( z = z, w_bar = wghts, exp2norm2 = sum(beta_moments_2) )
  M2
  M2Analytic( z = z, w_bar = wghts, method = "Bayesian" )
  Vz_ml$cov / (n + 1)
  mean( apply( bb_mat, 2, var ) )
  
  
  # Skewness of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^3 ) } ) )
  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  M3Exact( z = z, 
           w_bar = wghts, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2) )
  # M3BBFUN( x = z, M2 = M2, wghts = wghts )
  # M3BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, skewFUN ) )
  
  
  # Kurtosis of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^4 ) } ) )
  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  # M4BBFUN( x = z, M2 = M2, wghts = wghts )
  # M4BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, kurtFUN ) )  
  
  
 
  
  
  
  # --------------------------------------------------------------------------
  # Case 2b:  Symmetric alpha, but sum(alpha) < n 
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  alpha <- rep(1/n, n)  # --> leads to a scaling factor for the variance of 0.5
  wghts <- alpha / sum(alpha)
  w_bar <- wghts
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected squared 2-norm
  expected_2normSq <- mean( unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) ) )
  
  # Exact expected (2-norm)^2, (3-norm)^3, (4-norm)^4
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  Vz_ml$cov; sum( wghts * (z - sum(wghts*z))^2 )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
  # Variance of sample mean
  M2 <- M2Exact( z = z, w_bar = wghts, exp2norm2 = sum(beta_moments_2) )
  M2
  M2Analytic( z = z, w_bar = wghts, method = "Bayesian" )
  Vz_ml$cov / (n + 1)
  mean( apply( bb_mat, 2, var ) )
  
  
  # Skewness of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^3 ) } ) )
  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  M3Exact( z = z, 
           w_bar = wghts, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2) )
  # M3BBFUN( x = z, M2 = M2, wghts = wghts )
  # M3BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, skewFUN ) )
  
  
  # Kurtosis of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^4 ) } ) )
  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  # M4BBFUN( x = z, M2 = M2, wghts = wghts )
  # M4BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, kurtFUN ) )  
  
  
  
  
  

  
  # --------------------------------------------------------------------------
  # Case 3: Asymmetric alpha with sum(alpha) <> n 
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  alpha_n <- 1:n / sum(1:n) * n
  lambda <- 0.1
  alpha <- alpha_n * lambda
  wghts <- alpha / sum(alpha)
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  lNorm <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP[[j]], 1, lpNorm, p = 2 )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected squared 2-norm
  expected_2normSq <- mean( unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) ) )
  
  # Exact expected (2-norm)^2, (3-norm)^3, (4-norm)^4
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
 
  # Variance of sample mean
  M2 <- M2Exact( z = z, w_bar = wghts, exp2norm2 = sum(beta_moments_2) )
  M2
  M2Analytic( z = z, w_bar = wghts, method = "Bayesian" )
  Vz_ml$cov / (n + 1)
  mean( apply( bb_mat, 2, var ) )
  
  
  # Skewness of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^3 ) } ) )
  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  M3Exact( z = z, 
           w_bar = wghts, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2) )
  # M3BBFUN( x = z, M2 = M2, wghts = wghts )
  # M3BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, skewFUN ) )
  
  
  # Kurtosis of sample mean
  mean( apply( bb_mat, 2, FUN = function(x) { mean( (x - mean(x))^4 ) } ) )
  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  # M4BBFUN( x = z, M2 = M2, wghts = wghts )
  # M4BBFUN( x = z, M2 = NULL, wghts = wghts )
  # mean( apply( bb_mat, 2, kurtFUN ) )  
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 4a: FEV-Bias, symmetric
  # --------------------------------------------------------------------------
  
  
  
  
  # Dirichlet parameter
  alpha <- rep(2, n)
  wghts <- alpha / sum(alpha)
  
  # FEV-Bias
  fev_bias <- 20
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- lP_fev <- list()
  lNorm2 <- lNorm3 <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lP_fev[[j]] <- fevBias( x = lP[[j]], q = fev_bias )
    lNorm2[[j]] <- apply( lP_fev[[j]], 1, lpNorm, p = 2 )
    lNorm3[[j]] <- apply( lP_fev[[j]], 1, lpNorm, p = 3 )
    bb_mat[ ,j] <- apply( lP_fev[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP_fev, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected 2-norm
  expected_2normSq <- unlist( lapply( lNorm2, FUN = function(x) { mean(x^2) } ) )
  expected_3norm3 <- unlist( lapply( lNorm2, FUN = function(x) { mean(x^3) } ) )
  
  # Exact expected 2-norm
  # ... does not apply
  
  
  # Weighted sample variance
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   # wt = wghts, 
                   wt = w_mu,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                   # wt = wghts, 
                   wt = w_mu,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "unbiased" )
  
  # 2-norm of weights
  # lp_wbar <- lpNorm( x = wghts, p = 2 )
  # lp_wbar <- lpNorm( x = mu, p = 2 )
  lp_wbar <- apply( mumat, 2, lpNorm, p = 2 )
  
  # Compute variance of sample mean
  tmp <- cbind( as.numeric(Vz_ml$cov) * (1 + (expected_2normSq - 1) / (1 - lp_wbar^2) ),
                apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
  tmp
  apply( tmp, 2, mean )
  M2 <- mean(tmp[ ,1])
  
  
  
  
  V1 <- sum(w_mu)
  V2 <- sum(w_mu^2)
  V3 <- sum(w_mu^3)
  scl2 <- V1^2 / (V1^2 - V2)
  scl3 <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
  
  m2 <- sum( w_mu * (z - sum(w_mu * z))^2 )
  M2Exact( z = z, w_bar = w_mu, exp2norm2 = mean(expected_2normSq) )
  mean( apply( bb_mat, 2, momFUN, k = 2, scl = FALSE ) )  # fairly close
  
  
  m3 <- sum( w_mu * (z - sum(w_mu * z))^3 )
  m3 * scl3 * ( mean(expected_3norm3) - 3 * mean(expected_2normSq) * sum(w_mu^2) + 2 * sum(w_mu^3) )
  mean( apply( bb_mat, 2, momFUN, k = 3, scl = FALSE ) )  # fairly close
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 4b: FEV-Bias, Asymmetric
  # --------------------------------------------------------------------------
  
  # Dirichlet parameter
  alpha <- 1:n
  wghts <- alpha / sum(alpha)
  
  # FEV-Bias
  fev_bias <- 14
  
  # Run sampler
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- lP_fev <- list()
  lNorm2 <- lNorm3 <- list()
  lNorm2_rel <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lP_fev[[j]] <- fevBias( x = lP[[j]], q = fev_bias )
    lNorm2[[j]] <- apply( lP_fev[[j]], 1, lpNorm, p = 2 )
    lNorm3[[j]] <- apply( lP_fev[[j]], 1, lpNorm, p = 3 )
    w_mu <- apply( lP_fev[[j]], 2, mean )
    lNorm2_rel[[j]] <- apply( lP_fev[[j]], 1, function(x) { sum( (x - w_mu)^2 ) } )
    bb_mat[ ,j] <- apply( lP_fev[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Mean vector of the Dirichlet samples
  mumat <- do.call( cbind, lapply( lP_fev, FUN = function(P) { apply(P, 2, mean) } ) )
  w_mu <- apply( mumat, 1, mean )
  boxplot( t(mumat) )
  points( w_mu, col = 2 )
  points( wghts, col = 3 )
  
  
  # Expected 2-norm
  expected_2normSq <- unlist( lapply( lNorm2, FUN = function(x) { mean(x^2) } ) )
  expected_3norm3 <- unlist( lapply( lNorm3, FUN = function(x) { mean(x^3) } ) )
  expected_2normSq_rel <- mean( unlist( lNorm2_rel) )

    
  # Exact expected 2-norm
  # ... does not apply
 
  
  # Weighted sample variance
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   # wt = wghts, 
                   wt = w_mu,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                    # wt = wghts, 
                    wt = w_mu,
                    cor = FALSE, 
                    center = TRUE, 
                    method = "unbiased" )
  
  
  FUN <- function(i)
  {
    cov.wt( x = matrix(z, ncol = 1), 
            wt = mumat[ ,i],
            cor = FALSE, 
            center = TRUE, 
            method = "ML" )$cov
  }
  vz_ml <- unlist( lapply( 1:ncol(mumat), FUN = FUN ) )
  
  # 2-norm of weights
  # lp_wbar <- lpNorm( x = wghts, p = 2 )
  # lp_wbar <- lpNorm( x = mu, p = 2 )
  lp_wbar <- apply( mumat, 2, lpNorm, p = 2 )
  
  # Compute variance of sample mean
  tmp <- cbind( as.numeric(Vz_ml$cov) * (1 + (expected_2normSq - 1) / (1 - lp_wbar^2) ),
                apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
  tmp
  apply( tmp, 2, mean )
  M2 <- mean(tmp[ ,1])
  
  
  
  V1 <- sum(w_mu)
  V2 <- sum(w_mu^2)
  V3 <- sum(w_mu^3)
  scl2 <- V1^2 / (V1^2 - V2)
  scl3 <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
  

  Vz_unb$cov * expected_2normSq_rel 
  M2Exact( z = z, w_bar = w_mu, exp2norm2 = mean(expected_2normSq) )
  mean( apply( bb_mat, 2, momFUN, k = 2, scl = FALSE ) )  # fairly close
  
  
  
  
  m3 <- sum( w_mu * (z - sum(w_mu * z))^3 )
  m3 * scl3 * ( mean(expected_3norm3) - 3 * mean(expected_2normSq) * sum(w_mu^2) + 2 * sum(w_mu^3) )
  mean( apply( bb_mat, 2, momFUN, k = 3, scl = FALSE ) )  # fairly close
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
