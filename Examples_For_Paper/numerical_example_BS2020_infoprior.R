  
  
  ############################################################################
  ### NUMERICAL EXAMPLE OF BERTSIMAS AND STURT 2020 REVISITED 2
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.09.2022
  # First version:    16.06.2022
  # --------------------------------------------------------------------------
  
  # Example in Bertsimas and Sturt (2020) - informative prior
  
  
  # Structure:
  #
  # 1) Asymmetry
  # 2) Concentration
  # 3) Constrained domain - shadow Dirichlet - monotonic pmf's
  # 4) Mixture
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(orthopolynom)
  require(PDQutils)
  require(rgl)
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  # --------------------------------------------------------------------------
  # Load data
  # --------------------------------------------------------------------------
  
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
  boot_mat <- env$boot_mat
  bb_mat <- env$bb_mat
  
  z_demeaned <- scale(z, TRUE, FALSE)
  
  
  
  # --------------------------------------------------------------------------
  # 1) Asymmetry
  # --------------------------------------------------------------------------
  
  alpha <- 1:n / sum(1:n) * n
  # alpha <- rep(1, n)
  bb_mat_asy <- matrix(NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    samples <- rdirichlet( n = B, alpha = alpha ) 
    bb_mat_asy[ ,j] <- apply( samples, 1, function(w) { sum(w * z) } )
  }
  mr_bb_asy <- matrix(NA, nrow = n_boot, ncol = 20 )
  mr_bb_asy[ ,1] <- apply( bb_mat_asy, 2, mean ) 
  for ( i in 2:ncol(mr_bb_asy) ) {
    mr_bb_asy[ ,i] <- apply( bb_mat_asy, 2, FUN = function(x) { mean((x - mean(x))^i) } ) 
  } 

  
  
  # Exact moments
  w_bar <- alpha / sum(alpha)
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
  
  m1 <- sum( w_bar * z )
  m2 <- M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments_2) )
  m2b <- M2Analytic( z = z, w_bar = w_bar, method = "Bayesian" )
  m3 <- M3Exact( z = z,
                 w_bar = w_bar,
                 exp3norm3 = sum(beta_moments_3),
                 expw2 = beta_moments_2 )
  m3b <- M3BBFUN( x = z, M2 = 0, wghts = w_bar, scl = FALSE )
  m4 <- M4BBFUN( x = z, M2 = 0, wghts = w_bar, scl = FALSE )
  
  
  
  mean( mr_bb_asy[ ,1] ); m1
  mean( mr_bb_asy[ ,2] ); m2; m2b
  mean( mr_bb_asy[ ,3] ); m3; m3b
  mean( mr_bb_asy[ ,4] ); m4
  
    
  
  plot( density( mr_bb_asy[ ,1]) )
  abline( v = mean(mr_bb_asy[ ,1]) )
  
  plot( density( mr_bb_asy[ ,2]) )
  abline( v = mean( mr_bb_asy[ ,2]) )
  abline( v = m2, col = 2 )
  
  plot( density( mr_bb_asy[ ,3]) )
  abline( v = mean( mr_bb_asy[ ,3]) )
  abline( v = m3, col = 2 )
  abline( v = m3b, col = 3 )
  
  plot( density( mr_bb_asy[ ,4]) )
  abline( v = mean( mr_bb_asy[ ,4]) )
  abline( v = m4, col = 2 )  
  

  # Cornish-Fisher expansion:
  # Approximate the quantile using PDQutils
  # based on the first four central moments

  moments <- c(0, m2, m3b, m4)
  moments <- c(0, apply( mr_bb_asy[ ,-1], 2, mean) )[ c(1:4) ]
  
  raw_cumulants <- moment2cumulant(moments)
  z0_approx_bb_asy <- qapx_cf( p = p_th, 
                               raw.cumulants = raw_cumulants, 
                               support = c(-Inf, Inf), 
                               lower.tail = TRUE, 
                               log.p = FALSE )
  z0_approx_bb_asy <- z0_approx_bb_asy + m1

  
  z0_bb_asy <- apply( bb_mat_asy, 2, quantile, p_th )
  range(z0_bb_asy)
  mean(z0_bb_asy)
  quantile( as.numeric(bb_mat_asy), p_th )
  z0_approx_bb_asy
  
  
  
  plot( density(z0_bb_asy) )
  abline( v = mean(z0_bb_asy) )
  abline( v = z0_approx_bb_asy, col = 2 )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # 2) Concentration
  # --------------------------------------------------------------------------
  
  lambda <- 0.1
  alpha_conc <- alpha * lambda
  bb_mat_conc <- matrix(NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    samples <- rdirichlet( n = B, alpha = alpha_conc ) 
    bb_mat_conc[ ,j] <- apply( samples, 1, function(w) { sum(w * z) } )
  }
  theta_conc <- as.numeric(bb_mat_conc)
  mr_bb_conc <- matrix(NA, nrow = n_boot, ncol = 20 )
  mr_bb_conc[ ,1] <- apply( bb_mat_conc, 2, mean ) 
  for ( i in 2:ncol(mr_bb_conc) ) {
    mr_bb_conc[ ,i] <- apply( bb_mat_conc, 2, FUN = function(x) { mean((x - mean(x))^i) } ) 
  } 
  
  
  # Exact moments
  w_bar <- alpha_conc / sum(alpha_conc)
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha_conc), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha_conc, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha_conc, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha_conc, FUN = betaFUN, m = 4 ) )
  
  m1 <- sum( w_bar * z )
  m2 <- M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments_2) )
  m2b <- M2Analytic( z = z, w_bar = w_bar, method = "Bayesian" )
  m3 <- M3Exact( z = z,
                 w_bar = w_bar,
                 exp3norm3 = sum(beta_moments_3),
                 expw2 = beta_moments_2 )
  m3b <- M3BBFUN( x = z, M2 = 0, wghts = w_bar, scl = FALSE )
  m4 <- M4BBFUN( x = z, M2 = 0, wghts = w_bar, scl = FALSE )
  
  
  sum(beta_moments_2)
  mean( apply( samples, 1, function(x) { sum(x^2) } ) )
  
  
  mean( m1_bb_conc ); m1
  mean( m2_bb_conc ); m2; m2b     # m2b is not correct (requires a0 = n)
  mean( m3_bb_conc ); m3; m3b
  mean( m4_bb_conc ); m4
  
  
  
  # Cornish-Fisher expansion:
  # Approximate the quantile using PDQutils
  # based on the first four central moments
  
  moments <- c(0, m2, m3, m4)
  moments <- c(0, apply( mr_bb_conc[ ,-1], 2, mean) )[ c(1:4) ]
  raw_cumulants <- moment2cumulant(moments)
  z0_approx_bb_conc <- qapx_cf( p = p_th, 
                               raw.cumulants = raw_cumulants, 
                               support = c(-Inf, Inf), 
                               lower.tail = TRUE, 
                               log.p = FALSE )
  z0_approx_bb_conc <- z0_approx_bb_conc + m1

  
  z0_bb_conc <- apply( bb_mat_conc, 2, quantile, 0.025 )
  range(z0_bb_conc)
  mean(z0_bb_conc)
  quantile( as.numeric(bb_mat_conc), p_th )
  z0_approx_bb_conc

  
  
  
  
  plot( density(z0_bb_conc) )
  abline( v = mean(z0_bb_conc) )
  abline( v = z0_approx_bb_conc, col = 2 )
  
  
  
  
  
  
  
  
  lTmp <- list()
  for ( i in 2:20 ) {
    lTmp[[i-1]] <- apply( bb_mat_conc, 2, FUN = function(x) { sum((x - mean(x))^i) / length(x) } ) 
  }
  m1_bb_conc <- apply( bb_mat_conc, 2, mean ) 
  
  
  moments <- unlist( lapply( lTmp, FUN = mean ) ) 
  raw_cumulants <- moment2cumulant(c(0, moments))
  z0_approx_bb_conc <- qapx_cf( p = p_th, 
                                raw.cumulants = raw_cumulants, 
                                support = c(-Inf, Inf), 
                                lower.tail = TRUE, 
                                log.p = FALSE )
  z0_approx_bb_conc <- z0_approx_bb_conc + mean(m1_bb_conc)
  
  range(z0_bb_conc)
  mean(z0_bb_conc)
  z0_approx_bb_conc
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # 3) Constrained domain - shadow Dirichlet - monotonic pmf's
  # --------------------------------------------------------------------------
  
  M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
  
  alpha <- rep(1, n)
  # alpha <- 1:n / sum(1:n) * n
  w_bar <- alpha / sum(alpha)
  bb_mat_SD <- matrix(NA, nrow = B, ncol = n_boot)
  for ( j in 1:n_boot ) {
    samples <- rdirichlet( n = B, alpha = alpha )
    shadow <- t( apply( samples, 1, function(w) { M %*% w } ) ) 
    # bb_mat_SD[ ,j] <- apply( samples, 1, function(w) { sum(w * z) } )
    bb_mat_SD[ ,j] <- apply( shadow, 1, function(w) { sum(w * z) } )
  }
  
  mr_bb_SD <- list( apply( bb_mat_SD, 2, mean ) )
  for ( i in 2:20 ) {
    mr_bb_SD[[i]] <- apply( bb_mat_SD, 2, function(x) { mean((x - mean(x))^i) } )
  }
  
  moments_bb <- c(0, unlist( lapply( mr_bb_SD, mean ) ) )
  theta <- as.numeric(bb_mat_SD)
  moments_bb <- c(0, unlist( lapply( 2:20, FUN = function(i) { sum((theta - mean(theta))^i) / length(theta) } ) ))
  moments_bb
  moments_bb <- moments_bb[1:4]
  
  
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  beta_moments_4 <- unlist( lapply( alpha, FUN = betaFUN, m = 4 ) )
    
  z_shadow <- t(M) %*% z
  m1 <- sum(w_bar * z_shadow)
  m2 <- M2Exact( z = z, 
                 w_bar = w_bar,
                 M = M,
                 exp2norm2 = sum(beta_moments_2) )
  m3 <- M3Exact( z = z, 
                 w_bar = w_bar,
                 M = M,
                 exp3norm3 = sum(beta_moments_3),
                 expw2 = beta_moments_2 )
  
  m1; mean(theta)
  m2; mean(mr_bb_SD[[2]]); var(theta)
  var(z_shadow) * (sum(beta_moments_2) - sum(w_bar^2))
  m3; mean(mr_bb_SD[[3]])
  
  Vz_unb <- cov.wt( x = z_shadow,
                    wt = w_bar,
                    method = "unbiased" )
  Vz_unb$cov
  
  
  
  moments_bb <- c(0, unlist(lapply(mr_bb_SD, mean))[-1])
  moments_exact <- c(0, m2, m3)
  cbind( moments_exact, moments_bb[1:3] )
  
  moments <- moments_exact
  # moments <- c(moments_exact, moments_bb[-c(1:3)][1])
  # moments <- moments_bb[1:4]
  
  raw_cumulants <- moment2cumulant( moments )
  z0_approx_bb_SD <- qapx_cf( p = p_th, 
                              raw.cumulants = raw_cumulants, 
                              support = c(-Inf, Inf), 
                              lower.tail = TRUE, 
                              log.p = FALSE )
  z0_approx_bb_SD <- z0_approx_bb_SD + mean( mr_bb_SD[[1]] )
  
  
  z0_bb_SD <- apply( bb_mat_SD, 2, quantile, p_th )
  range(z0_bb_SD)
  mean(z0_bb_SD); quantile( theta, p_th )
  z0_approx_bb_SD
  
  
  
  # -----
  # Varsi
  # -----
  
  # q_vec <- seq( from = 1e-12, to = 0.999, length.out = 10^4 )
  # z_quantile <- quantile(z, q_vec)
  # z_shadow_quantile <- quantile(z_shadow, q_vec)
  z_shadow_quantile <- seq( from = mean(z_shadow)-sd(z_shadow), 
                            to = mean(z_shadow)+sd(z_shadow), length.out = 10^4 )
  FUN <- function(z0) { tail(as.numeric(varsi( mu = z_shadow, b = z0)), 1) }
  p_vec_shadow <- unlist( lapply( z_shadow_quantile, FUN ) )
  
  plot( x = z_shadow_quantile, y = p_vec_shadow )

  
  quantile(theta, p_th)
  z_shadow_quantile[ which(p_vec_shadow > p_th)[1] ]
  
  
  # Find the quantile H(0.025) iteratively
  varsiFindQuantile( z = z_shadow, th = p_th )
  
  
  
  
  # Derivative
  idx <- 3:length(p_vec_shadow)
  d_vec <- (p_vec_shadow[idx] - p_vec_shadow[idx-2]) / abs(z_shadow_quantile[idx] - z_shadow_quantile[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  z_vec <- z_shadow_quantile[-c(1, length(z_shadow_quantile))]
  plot( x = z_vec, y = d_vec_std )
  

  # Variance from density
  sum( d_vec_std * (z_vec - sum(z_vec * d_vec_std))^2 )
  var(theta)
  M2Exact( z = z, M = M, exp2norm2 = 2/(n+1) )
  
  
  
  
  
  
  
  
  
  
  ######################
  m1 <- mean(z)
  m2 <- sum( (z - m1)^2 ) / (n * (n + 1))
  tmp <- unlist( lapply( 3:20, FUN = function(k) { attr( getMk( z = z, k = k ), "B" ) } ) )
  moments_getmk <- c(0, m2, tmp)
  cbind( moments_bb, moments_getmk )
    
  raw_cumulants <- moment2cumulant(moments_bb)
  raw_cumulants <- moment2cumulant(moments_getmk)
  z0_approx_bb_SD <- qapx_cf( p = p_th, 
                              raw.cumulants = raw_cumulants, 
                              support = c(-Inf, Inf), 
                              lower.tail = TRUE, 
                              log.p = FALSE )
  z0_approx_bb_SD <- z0_approx_bb_SD + m1
  
  
  z0_bb_SD <- apply( bb_mat_SD, 2, quantile, p_th )
  range(z0_bb_SD)
  mean(z0_bb_SD); quantile( theta, p_th )
  z0_approx_bb_SD
  
    