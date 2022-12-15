  
  
  ############################################################################
  ### EXACT MOMENTS - CLASSICAL BOOTSTRAP - WIP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.05.2021
  # First version:    22.05.2021
  # --------------------------------------------------------------------------
  
  
  # Show that the standard error of the ordinary bootstrap can be computed 
  # exactly.
  
  
  
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
  
  
 
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 10^2
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1111)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  
  # x <- rnorm( n = n )
  # z <- x^2
  
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  
  tic <- Sys.time()
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
  
  
  
  
  # # --------------------------------------------------------------------------
  # # Varsi, i.e., exact BB
  # # --------------------------------------------------------------------------
  # 
  # se <- mean( apply(boot_mat, 2, sd) )
  # 
  # 
  # z_vec <- seq( from = 0, to = mean(z) + se * 10, length.out = 10^5 )
  # # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  # FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) }
  # tic <- Sys.time()
  # p_vec <- try( unlist( lapply( z_vec, FUN ) ) ) 
  # if ( inherits(p_vec, "try-error") ) {
  #   q_vec <- seq( from = 0.1, to = 0.999, length.out = 10^5 )
  #   z_vec <- quantile(z, q_vec)
  #   p_vec <- try( unlist( lapply( z_vec, FUN ) ) ) 
  # }
  # (toc_varsi <- Sys.time() - tic)
  # 
  # # Derivative
  # idx <- 3:length(p_vec)
  # # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(z_vec[idx] - z_vec[idx-2])
  # d_vec_std <- d_vec / sum(d_vec)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Standard error
  # --------------------------------------------------------------------------
  
  # Ordinary bootstrap
  stderr_boot <- apply( boot_mat, 2, sd )
 
  stderr_exact_boot <- sqrt( var(z) * (n-1) / n / n )
  mean(stderr_boot)
  stderr_exact_boot

  
  # Bayesian bootstrap (BB)
  stderr_bb <- apply( bb_mat, 2, sd )
 
  # Exact BB
  stderr_exact_bb <- sqrt( var(z) * ( (2 / (n + 1)) - 1 / n ) )
  stderr_exact_cales <- sqrt( sum((z - mean(z))^2) / (n * (n + 1)) )  # same same
 
  # # Exact BB by Varsi
  # z_eval <- z_vec[-c(1, length(z_vec))]
  # sum( d_vec_std * z_eval ); mean(z)
  # stderr_exact_bb_varsi <- sqrt( sum( d_vec_std * (z_eval)^2 ) - sum( d_vec_std * z_eval )^2 )
  
  
  
  # Compare
  SE <- setNames( c(mean(stderr_boot), 
                    mean(stderr_bb), 
                    stderr_exact_bb,
                    # stderr_exact_bb_varsi,
                    sd(z) / sqrt(n)),
                  c("avg_boot",
                    "avg_bb", 
                    "exact_bb",
                    # "exact_bb_varsi",
                    "normal") )
  SE
  
  
  
  
  # --------------------------------------------------------------------------
  # Skewness
  # --------------------------------------------------------------------------
 
  # Ordinary bootstrap
  M3_boot <- apply( boot_mat, 2, skewFUN )
  
  sk_biased <- skewFUN( x = z )
  M3_exact_boot <- sk_biased * n / n^(3/2)
  sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2)
  mean(M3_boot)
  
  # Bayesian bootstrap (BB)
  M3_bb <- apply( bb_mat, 2, FUN = skewFUN )
  
  # Exact BB
  M3_exact_bb <- M3FUN( x = z, M2 = stderr_exact_bb^2 )
  
  
  # Compare
  M3 <- setNames( c(mean(M3_boot), 
                    M3_exact_boot,
                    mean(M3_bb), 
                    M3_exact_bb,
                    0),
                  c("avg_boot",
                    "exact_boot",
                    "avg_bb", 
                    "exact_bb",
                    "normal") )
  M3
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Kurtosis
  # --------------------------------------------------------------------------
  
  # Ordinary bootstrap
  M4_boot <- apply( boot_mat, 2, kurtFUN )

  # Exact ordinary bootstrap
  kurt_biased <- kurtFUN( x = z )
  M4_exact_boot <- kurt_biased * n / n^2
  sum( (z - mean(z))^4 ) / sum( (z - mean(z))^2 )^2
  
  # Bayesian bootstrap (BB)
  M4_bb <- apply( bb_mat, 2, FUN = kurtFUN )
  
  # Exact BB
  M4_exact_bb <- M4FUN( x = z, M2 = stderr_exact_bb^2 )
  
  
  # Compare
  M4 <- setNames( c(mean(M4_boot), 
                    M4_exact_boot,
                    mean(M4_bb), 
                    M4_exact_bb,
                    # kurtFUN( rnorm(10^7, mean(z), sd(z) / sqrt(n)))),
                    3),
                  c("avg_boot",
                    "exact_boot",
                    "avg_bb", 
                    "exact_bb",
                    "normal") )
  M4
  
  
  SE
  M3
  M4
  
  
  
  df <- 10
  x_tmp <- rt( n = 10^7, df = df )
  kurtFUN( x_tmp )
  (3 * df - 6) / (df - 4) # Analytic kurtosis of t-dist
  
  plot( density(x_tmp) ) 
  lines( density(rnorm(10^7)), col = 2)
  
  
  
  
  # --------------------------------------------------------------------------
  # n -> inf
  # --------------------------------------------------------------------------
  
  # Notes: 
  # Finite n, we have M2 < sd(z) / sqrt(n), but non-zero skew and excess kurtosis
  # As n --> inf, skew -> 0 and kurt -> 3 (both, z \sim Normal and z \sim student-t)
  
  n <- 10^2
  z <- rnorm(n)
  M2 <- sqrt( var(z) * ( (2 / (n + 1)) - 1 / n ) )
  M3FUN( x = z, M2 = M2^2 )
  M4FUN( x = z, M2 = M2^2 )
  
  M2; sd(z) / sqrt(n)
  
  
  n <- 10^7
  df <- 10
  z <- rt( n = n, df = df )
  M2 <- sqrt( var(z) * ( (2 / (n + 1)) - 1 / n ) )
  M3FUN( x = z, M2 = M2^2 )
  M4FUN( x = z, M2 = M2^2 )
  (3 * df - 6) / (df - 4) 
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact standard error identity
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(gtools)
  
  
  n <- 3
  seed <- 1234
  
  P(n)
  comb( 2*n-1, n )
  
  # Random variable
  set.seed(seed)
  x <- rt( n = n, df = 5 )
  z <- x^3
  Vz <- var(z)
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), wt = rep(1/n, n), center = TRUE, method = "unbiased" )$cov
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), wt = rep(1/n, n), center = TRUE, method = "ML" )$cov
  Vz; Vz_unb; Vz_ml
  
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # 2-norm
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  
  
  # permutations
  lP <- list()
  lhs_rhs_same_same <- rep( NA, ncol(parts_mat) )
  lhs_lhs2_same_same <- rep( NA, ncol(parts_mat) )
  for( j in 1:ncol(parts_mat) ) {
    
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    lP[[j]] <- perm_unique
    theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
    LHS <- (lp_parts[j]^2 - 1/n) * Vz
    LHS2 <- lp_parts[j]^2 * Vz_ml
    lp_parts[j]^2 * Vz - Vz / n
    lp_parts[j]^2 * Vz - (Vz - Vz_ml)
    #
    Vz * (lp_parts[j]^2 - 1) + Vz_ml  # !!!!!
    Vz_ml * (1 + (lp_parts[j]^2 - 1) / (1 - 1/n))   # !!!!!!!!
    Vz_unb * (1 - 1/n) == Vz_ml
    #
    RHS <- sum( (theta - mean(theta))^2 ) / length(theta)
    lhs_rhs_same_same[j] <- abs( LHS - RHS ) < 1e-12
    lhs_lhs2_same_same[j] <- abs( LHS - LHS2 ) < 1e-12
  }
  lhs_rhs_same_same
  lhs_lhs2_same_same
  
  
  
  Vz_ml * (1 + (lp_parts[j]^2 - 1) / (1 - 1/n))
  Vz_unb * (lp_parts[j]^2 - 1/n)  # same same
  
  
  
  k <- 2
  P <- lP[[k]]
  theta <- apply( P / n, 1, function(p) { sum(p * z) } )
  # LHS <- (lp_parts[j]^2 - 1/n) * Vz
  LHS <- Vz_ml * ( 1 + (lp_parts[k]^2 - 1) / (1 - 1/n) )
  RHS <- sum( (theta - mean(z))^2 ) / length(theta)
  LHS; RHS
  
  theta
  a <- sum( (z - mean(z))^2 )
  b <- sum( (theta - mean(z))^2 )
  a_scl <- (1 + (lp_parts[k]^2 - 1) / (1 - 1/n)) / n
  b_scl <- length(theta)
  a
  b
  a * a_scl
  b * (1 / b_scl)
  mean(theta^2) - mean(z)^2
  (mean(z^2) - mean(z)^2) * (1 + (lp_parts[k]^2 - 1) / (1 - 1/n))
  f <- (lp_parts[k]^2 - 1) / (1 - 1/n)
  mean(z^2) + mean(z^2) * lp_parts[k]^2 / (1 - 1/n) - mean(z^2) / (1 - 1/n) - mean(z)^2 - mean(z)^2 * f
  
  
  mean(theta^2)
  (mean(z^2) - mean(z)^2) * (lp_parts[k]^2 - 1) / (1 - 1/n) + mean(z^2)
  
  
  B <- t( apply( P, 1, function(p) { p * z } ) )
  B
  
  rho <- parts_mat[ ,k] / n
  c(rho^2, -1) / (1 - 1/n) 
  sum(c(rho^2, -1)) / (1 - 1/n) 
  sum( c(rho^2, -1) / (1 - 1/n) )
  (sum(rho^2) - 1) / (1  - 1/n)
  
  A <- matrix(z^2, ncol = 1) %*% matrix( c(rho^2, -1) / (1 - 1/n) * 1/n, nrow = 1)
  
  
  
  lp_parts[k]^2
  sum(rho^2)
  sum((rho * n)^2) / (n * nrow(P))
  sum( apply( P, 1, function(p) { sum(p * z) } ) ) / (nrow(P) * n)
  sum( apply( P / n, 1, function(p) { sum(p * z) } ) ) / (nrow(P) * n)
  mean(theta) / n
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact standard error identity  - informative prior
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(gtools)
  
  
  n <- 8
  seed <- 1234
  n_boot <- 50
  B <- 10^4
  
  P(n)
  comb( 2*n-1, n )
  
  # Weights
  wghts <- 1:n / sum(1:n)

  # Random variable
  set.seed(seed)
  x <- rt( n = n, df = 5 )
  z <- x^3
  mu <- sum( wghts * z )
  s <- var(z)
  s_unb <- cov.wt( x = matrix(z, ncol = 1), wt = wghts, center = TRUE, method = "unbiased" )$cov
  s_ml <- cov.wt( x = matrix(z, ncol = 1), wt = wghts, center = TRUE, method = "ML" )$cov
  s; s_unb; s_ml
  
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # 2-norm
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  lp_wghts <- lpNorm( wghts, p = 2 )
  
  # permutations
  lP <- lprob <- lTheta <- list()
  lhs_rhs_same_same <- rep( NA, ncol(parts_mat) )
  lhs <- lhs_rhs_same_same
  rhs <- lhs
  for( j in 1:ncol(parts_mat) ) {
    
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    prob_perm <- apply( perm_unique, 1, dmultinom, prob = wghts, log = FALSE )
    theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
    LHS <- (lp_parts[j]^2 - lp_wghts^2) * s_unb
    s_unb * (lp_parts[j]^2 - 1) + s_ml  # !!!!!
    s_ml * (1 + (lp_parts[j]^2 - 1) / (1 - lp_wghts^2))   # !!!!!!!!
    s_unb * (1 - lp_wghts^2) == s_ml
    RHS <- sum( prob_perm / sum(prob_perm) * (theta - mu)^2 )
    # RHS <- sum( prob_perm * (theta - mu)^2 )
    # cov.wt( matrix(theta, ncol = 1), 
    #         wt = prob_perm, 
    #         center = mu, 
    #         method = "ML")$cov 
    lhs_rhs_same_same[j] <- abs( LHS - RHS ) < 1e-12
    lP[[j]] <- perm_unique
    lprob[[j]] <- prob_perm
    lTheta[[j]] <- theta
    lhs[j] <- LHS
    rhs[j] <- RHS
  }
  lhs_rhs_same_same
  
  mu
  sum( unlist( lapply( 1:length(lTheta), FUN = function(i) { sum( lTheta[[i]] * lprob[[i]] ) } ) ) )
  
  sum(rhs)
  
  
  
  
  # Expected 2-norm under prior probabilities
  lp_lP <- unlist( lapply( lP, function(P) { apply( P / n, 1, lpNorm, p = 2 ) } ) )
  elp2 <- sum( lp_lP^2 * unlist(lprob) )
  elp2
  
  
  # Compare with weighted bootstrap
  boot_mat_w <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        weights = wghts,
                        R = B )
    boot_mat_w[ ,j] <- Boot$t
  }
  
  mean( apply( boot_mat_w, 2, var ) )
  s_unb * ( elp2 - lp_wghts^2 )  # close enough
  sum( lhs * unlist( lapply( lprob, sum ) ) ) # same same 
  cov.wt( matrix(unlist(lTheta), ncol = 1), 
          wt = unlist(lprob), 
          center = TRUE, 
          method = "ML")$cov  # same same
 
  
  M2 <- s_ml * (1 + (elp2 - 1) / (1 - lp_wghts^2))
  A <- (1 + (elp2 - 1) / (1 - lp_wghts^2))
  A
  M2 / s_ml
  1 / A
  elp2
  2/n - 1/n^2
  M2
  s_ml / n
  
  
  
  
  
  # Question:
  # Is E(||w||_2^2) - ||\bar{w}||_2^2 = E(||w - \bar{w}||_2^2) constant for 
  # fixed n and for any distribution of w \sim \Delta?
  # No. 
  # But it is for asymmetric Dirichlet if sum(alpha) = n.
  
  # Does it hold when using FEV-Bias?
  
  
  
  n <- 10
  seed <- 1234
  n_boot <- 50
  B <- 10^5
  df <- 4
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  m <- 1; wghts <- (1:n)^m / sum((1:n)^m)
  sum(wghts)
  # wghts <- rep(1/n, n)
  alpha <- wghts * n
  # alpha <- wghts * 0.05
  # alpha <- wghts * n * 10000
  
  # Bayesian bootstrap
  
  P <- rdirichlet( n = B, alpha = alpha )
  lp2_w <- apply( P, 1, lpNorm, p = 2 )
  lp2_wbar <- lpNorm( x = wghts, p = 2 )
  lp2_wminuswbar <- apply( P, 1, lpNorm, p = 2, x_bar = wghts )
  
  mean(lp2_w^2) - lp2_wbar^2
  mean(lp2_wminuswbar^2)
  
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  sum(beta_moments) - lpNorm( x = wghts, p = 2 )^2
  
  theta <- apply( P, 1, function(p) { sum(p * z) } )
  mean(theta)
  var(theta)
  
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
  

  sum( (theta - mean(theta))^2 ) / (length(theta) - 1)
  s_unb_w * (sum(beta_moments) - lpNorm( x = wghts, p = 2 )^2)
  s_unb_w * (mean(lp2_w^2) - lp2_wbar^2)
  s_unb_w * mean(lp2_wminuswbar^2)
  s_ml_w * (n + 1)^-1
 
  (1 + (2/(n+1) - 1) / (1 - 1/n))
  (n + 1)^-1
  (1 + (sum(beta_moments) - 1) / (1 - lpNorm( x = wghts, p = 2 )^2))
  ( sum(beta_moments) - lpNorm( x = wghts, p = 2 )^2 ) / (1 - lpNorm( x = wghts, p = 2 )^2)
  
  
  
  
  
  # Ordinary bootstrap
  
  
  
  
  
  
  
  
  
  
  
  

  # --------------------------------------------------------------------------
  # Exact skewness and kurtosis identity
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
  w_bar <- rep(1 / n, n)
  P <- rdirichlet( n = B, alpha = rep(1, n) )
  P_lp2 <- apply( P, 1, lpNorm, p = 2 )
  P_lp3 <- apply( P, 1, lpNorm, p = 3 )
  
 
  # Ordinary bootstrap
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  tic <- Sys.time()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  var_boot <- apply( boot_mat, 2, var )
  sk_boot <- apply( boot_mat, 2, skewFUN )
  kurt_boot <- apply( boot_mat, 2, kurtFUN )
  cm3_boot <- apply( boot_mat, 2, function(x) { sum( (x - mean(x))^3 ) / length(x) } )
  
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # Reorder results by decreasing total probability
  ordering <- rev( order( prob_tot) ) 
  parts_mat <- parts_mat[ ,ordering]
  draws <- draws[ ,ordering]
  prob <- prob[ordering]
  prob_tot <- prob_tot[ordering]
  
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  
  
  # permutations
  lP <- list()
  lTheta <- list()
  lLHS <- lRHS <- list()
  lLHS3 <- lRHS3 <- list()
  lRHS4 <- list()
  
  for( j in 1:ncol(parts_mat) ) {
    
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    lP[[j]] <- perm_unique
    theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
    lTheta[[j]] <- theta
    lLHS[[j]] <- (lp_parts[j]^2 - 1/n) * s_unb
    lRHS[[j]] <- sum( (theta - mean(theta))^2 ) / length(theta)
    
    lLHS3[[j]] <- (lp_parts[j]^2 - 1/n) * s_unb
    
    lRHS3[[j]] <- sum( (theta - mean(theta))^3 ) / length(theta)
    lRHS4[[j]] <- sum( (theta - mean(theta))^4 ) / length(theta)
  }

  # Weights for each theta
  lW <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  w <- unlist(lW)
  sum(w)
  
  plot( unlist(lRHS), unlist(lLHS) )
  plot( unlist(lRHS3), unlist(lLHS3) )
  

  
  
  
  
  # Mean
  sum( w * unlist(lTheta) )
  mean(z) # same same
  unlist( lapply( lTheta, FUN = mean ) )
  
  M1 <- mean(z)
  
  
  # Variance
  cov.wt( x = matrix(unlist(lTheta), ncol = 1), 
          wt = w, 
          cor = FALSE, 
          center = TRUE, 
          method = "ML" )$cov 
  sum( w * ( unlist(lTheta) - sum( w * unlist(lTheta) ) )^2 ) # same same
  s_unb * ( sum(lp_parts^2 * prob_tot) - 1/n ) # same same
  s_ml * ( 1 + (sum(lp_parts^2 * prob_tot) - 1) / (1 - 1/n) ) # same same
  sum( unlist(lRHS) * prob_tot ) # same same
  sum( unlist(lLHS) * prob_tot ) # same same
  s_ml / n # same same
  mean(var_boot) # similar
  
  M2 <- s_ml / n
  
  plot( x = unlist(lRHS), y = lp_parts^2 )
  
  
  
  # Skewness
  
  skew.wt( x = unlist(lTheta),
           wt = w )
  sum( (z - mean(z))^3 ) / (n^3 * M2^(3/2) ) # same same
  sk * (n / n^(3/2)) # same same
  sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2) # same same
  mean(sk_boot) # not too different
  

  M3 <- sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2)
  
  
  A <- sum( (z - mean(z))^3 ) / (M3 * M2^(3/2))
  sum( (z - mean(z))^3 ) / (A * M2^(3/2))
  A
  n^3
  
  
  
  lRHS3
  tmp <- unlist( lapply( lRHS3, FUN = function(x) { x / M2^(3/2) } ) )
  sum( tmp * prob_tot )
  
  
  a <- sum( (z - mean(z))^3 ) / n^3
  a / M2^(3/2)
  sk * ( (1 / n) * (n / sqrt(n)) ) # same same
  
  
  
  
  # ---------------
  # Paper Angelova
  # See also Harald Cramer, Mathematical Methods of Statistics, p. 345
  # ---------------
  
  sum( (z - mean(z))^3 ) / n / n^2
  sum( (z - mean(z))^3 ) / n^3  # same same
  sum( prob_tot * unlist( lRHS3 ) )  # same same
  
  sum( (z - mean(z))^3 ) / n / n^2 / ( sum( (z - mean(z))^2 / n / n ) )^(3/2) # same same as M3
  M3  

  
  
  
 
  
  
  
  # Kurtosis
  
  kurt.wt( x = unlist(lTheta),
           wt = w )
  tmp <- unlist( lapply( lRHS4, FUN = function(x) { x / M2^(4/2) } ) )
  sum( tmp * prob_tot )
  
  M4 <- sum( tmp * prob_tot )

  
  A <- sum( (z - mean(z))^4 ) / (M4 * M2^2)
  sum( (z - mean(z))^4 ) / (A * M2^2)
  A
  
  
  AA <- M4 / (2 * ks + n)
  (2 * ks + n) * AA
  AA
  
  
  
  
  # ---------------
  # Paper Angelova
  # See also Harald Cramer, Mathematical Methods of Statistics, p. 345
  # ---------------

  M2 <- sum( (z - mean(z))^2 / n^2 )
  mu2 <- sum( (z - mean(z))^2 ) / n
  mu4 <- sum( (z - mean(z))^4 ) / n
  3 * mu2^2 / n^2 + (mu4 - 3 * mu2^2) / n^3
  3 * M2^2 - 3 * M2^2 / n + mu4 / n^3  # same same
  M2^2 * (3 - 3/n) + mu4 / n^3  # same same
  sum( prob_tot * unlist( lRHS4 ) )  # same same
  
  
  (3 * mu2^2 / n^2 + mu4 / n^3 - 3 * mu2^2 / n^3) / M2^2
  
  
  
  (3 * mu2^2 / n^2 + (mu4 - 3 * mu2^2) / n^3) / M2^2 # same same as M4
  (M2^2 * (3 - 3/n) + mu4 / n^3) / M2^2 # same same as M4
  3 - 3/n + mu4 / (n^3 * M2^2) # same same as M4
  M4
  
  
  sum( (z - mean(z))^4 ) / (mu2^2 * n^3) + 3 - 3 / n
  mu4 / (mu2^2 * n^2) + 3 - 3 / n
  
  
 
  
  
  
  
  
  
  
  
  
  
  # lTheta <- lapply( lP, FUN = function(P) { apply( P / n, 1, function(p) { sum(p * z) } ) } )
  cbind( unlist(lLHS), unlist(lRHS) )
  sum( unlist(lRHS) * prob_tot )
  s * ( sum(lp_parts^2 * prob_tot) - 1/n )  # same same
  
  rhs3 <- unlist(lRHS3)
  rhs3[is.na(rhs3)] <- 0
  sum( rhs3 * prob_tot )
  
  tmp <- lapply( 1:length(lTheta), FUN = function(i) { lRHS3[[i]] / (length(lTheta[[i]]) * lRHS[[i]]^(3/2)) })
  tmp <- unlist(tmp)
  tmp[is.na(tmp)] <- 0
  sum( tmp * prob_tot )
  
  
  
 
  
    
  
  
  # --------
  # Skewness
  # --------
  
  
  lRHS3
  rhs3 <- unlist(lRHS3)
  rhs3[is.na(rhs3)] <- 0
  sum( rhs3 * prob_tot )
  mean(cm3_boot)
  mean(cm3_boot) / mean(var_boot)^(3/2)
  
  
  
  sum( w * (unlist(lTheta) - mean(z))^3 ) / M2^(3/2)
  rhs3 <- unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(theta))^3 ) / length(theta) } ) )
  rhs3
  sum( rhs3 * prob_tot ) / M2^(3/2)  # same same
  sum( rhs3 / M2^(3/2) * prob_tot )  # same same
  sum( rhs3 * prob_tot ) / sum( prob_tot * unlist(lRHS) )^(3/2) # same same
  sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2) # same same 
  
  
  plot( x = rhs3 / M2^(3/2), y = (apply( parts_mat/n, 2, lpNorm, p = 3 )^3 - lpNorm(w_bar, p = 3)^3) / (lp_parts^2 - 1/n)^(3/2) )
  plot( x = rhs3 / M2^(3/2), y = apply( parts_mat/n, 2, lpNorm, p = 3 )^3 / (lp_parts^2)^(3/2) )
  
  
  
  sk <- skewFUN( x = z )
  sk
  sum( (z - mean(z))^3 ) / ( n * (sum( (z - mean(z))^2 ) / n)^(3/2) )
  
  ###
  sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2)
  ###
  
  A <- as.numeric(M3 / sk)
  A * sk; M3;  sum( (z - mean(z))^3 ) / ( n * (sum( (z - mean(z))^2 ) / n)^(3/2) ) * A
  
  A == n / n^(3/2) # different at machine error
  A - n / n^(3/2)
  ###
  A
  n / n^(3/2)
  ###
  
  
  
  
  # --------
  # Kurtosis
  # --------
  
  M4 <- kurt.wt( x = unlist(lTheta),
                 wt = w )
  krtsis <- kurtFUN( x = z )
  
  M4
  mean(kurt_boot) 

  krtsis
  sum( (z - mean(z))^4 ) / ( n * (sum( (z - mean(z))^2 ) / n)^(4/2) )
  
  A <- M4 / krtsis
  A
  n / n^2
  
  
  M4
  sum( (z - mean(z))^4 ) / (n^4 * M2^(4/2) )
  
  B <- sum( (z - mean(z))^4 ) / (M4 * M2^2)
  B
  n^4
  n / n^2
  
  
  
  
  
  
  
  gtools::permutations
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact skewness and kurtosis identity - loop
  # --------------------------------------------------------------------------
  
  
  B <- 10^4 + 1
  n_boot <- 50
  df <- 5
  n_vec <- 3:9
  
  m2mat <- matrix( NA, nrow = length(n_vec), ncol = 5,
                   dimnames = list(n_vec, c("wt", "exact", "boot", "bb", "s_biased")) )
  m3mat <- matrix( NA, nrow = length(n_vec), ncol = 7,
                   dimnames = list(n_vec, c("wt", "exact", "boot", "bb", "sk", "A", "B")) )
  m4mat <- matrix( NA, nrow = length(n_vec), ncol = 6,
                   dimnames = list(n_vec, c("wt", "boot", "bb", "ks", "A", "B")) )
  lZ <- list()
  
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    
    # Generate random sample
    set.seed(1234)
    x <- rt( n = n, df = df )
    z <- x^2
    lZ[[i]] <- z
    
    # Ordinary bootstrap
    boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
    tic <- Sys.time()
    for ( j in 1:n_boot ) {
      Boot <- boot::boot( data = z,
                          statistic = meanBoot,
                          R = B )
      boot_mat[ ,j] <- Boot$t
    }
    
    # Bayesian bootstrap
    bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
    lP <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = B, alpha = rep(1, n) )
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    # Compute sample moments on ordinary bootstrap means
    var_boot <- apply( boot_mat, 2, var )
    sk_boot <- apply( boot_mat, 2, skewFUN )
    kurt_boot <- apply( boot_mat, 2, kurtFUN )
    
    # Compute sample moments on bayesian bootstrap means
    var_bb <- apply( bb_mat, 2, var )
    sk_bb <- apply( bb_mat, 2, skewFUN )
    kurt_bb <- apply( bb_mat, 2, kurtFUN )
    
    # Partitions
    parts_mat <- partitions::parts( n )
    draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
    prob <- draws2prob( draws = draws )
    prob_tot <- prob * draws[2, ]
    
    # Reorder results by decreasing total probability
    ordering <- rev( order( prob_tot) ) 
    parts_mat <- parts_mat[ ,ordering]
    draws <- draws[ ,ordering]
    prob <- prob[ordering]
    prob_tot <- prob_tot[ordering]
    
    lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
    
    # permutations
    lP <- list()
    lTheta <- list()
    s <- var(z)
    s_biased <- sum( ( z - mean(z))^2 ) / n
    
    for( j in 1:ncol(parts_mat) ) {
      
      if ( length(unique(parts_mat[ ,j])) > 1 ) {
        perm <- gtools::permutations( n = n, 
                                      r = n, 
                                      v = parts_mat[ ,j], 
                                      set = FALSE,
                                      repeats.allowed = FALSE )
        perm_unique <- perm[!duplicated(perm), ]
      } else {
        perm_unique <- t(parts_mat[ ,j])
      }
      lP[[j]] <- perm_unique
      theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
      lTheta[[j]] <- theta
    }
    
    
    # Weights
    # use node activation probabilities as weights (for each node) and compute 
    # weighted moments on theta's.
    
    lW <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
    w <- unlist(lW)
    
    # Variance
    M2_wt <- cov.wt( x = matrix(unlist(lTheta), ncol = 1),
                     wt = w, 
                     cor = FALSE, 
                     center = TRUE, 
                     method = "ML")$cov
    M2_exact <- s * ( sum(lp_parts^2 * prob_tot) - 1/n )
    
    # Skewness
    sk <- skewFUN( x = z )
    M3_wt <- skew.wt( x = unlist(lTheta),
                      wt = w )
    M3_exact <- sum( (z - mean(z))^3 ) / (n^3 * M2_exact^(3/2) )
    
    # Kurtosis
    ks <- kurtFUN( x = z )
    M4_wt <- kurt.wt( x = unlist(lTheta),
                      wt = w )                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    
        
    # Store results
    m2mat[i, ] <- c(M2_wt,
                    M2_exact, 
                    mean(var_boot), 
                    mean(var_bb),
                    s_biased)
    m3mat[i, ] <- c(M3_wt, 
                    M3_exact, 
                    mean(sk_boot), 
                    mean(sk_bb),
                    sk,
                    sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2), 
                    n / n^(3/2) )
    m4mat[i, ] <-c(M4_wt, 
                   mean(kurt_boot), 
                   mean(kurt_bb),
                   ks,
                   sum( (z - mean(z))^4 ) / sum( (z - mean(z))^2 )^(4/2),
                   n / n^2)
    
    
  }
  
  m2mat
  cbind( m3mat, m3mat[ ,"sk"] * m3mat[ ,"B"] )
  cbind( m4mat, m4mat[ ,"ks"] * m4mat[ ,"B"], m4mat[ ,"wt"] / m4mat[ ,"ks"] )
  
  sum( (z - mean(z))^3 ) / ( m3mat[ ,"wt"] * m2mat[ ,"wt"]^(3/2) )
  n_vec^3
  
  sum( (z - mean(z))^4 ) / (m4mat[ ,"wt"] * m2mat[ ,"wt"]^2)
  
  
  FUN <- function(i) 
  {
    # sum( (z[1:(i+2)] - mean(z[1:(i+2)]))^3 ) / ( m3mat[i, "wt"] * m2mat[i, "wt"]^(3/2) )
    # sum( (lZ[[i]] - mean(lZ[[i]]))^3 ) / ( m3mat[i, "wt"] * m2mat[i, "wt"]^(3/2) )
    sum( (z[1:(i+2)] - mean(z[1:(i+2)]))^4 ) / ( m4mat[i, "wt"] * m2mat[i, "wt"]^(4/2) )
  }
  tmp <- unlist( lapply( 1:length(n_vec), FUN = FUN ) )
  tmp
  n_vec^3
  
  tmp
  
  
  
  m3_bb <- numeric(nrow(m2mat))
  m4_bb <- numeric(nrow(m2mat))
  for ( i in 1:nrow(m2mat) ) {
    m3_bb[i] <- M3FUN( x = lZ[[i]], M2 = m2mat[i, "wt"] )
    m4_bb[i] <- M4FUN( x = lZ[[i]], M2 = m2mat[i, "wt"] )
  }  
  m3_bb
  m4_bb
  
  cbind( m3mat, m3_bb )
  cbind( m4mat, m4_bb )
  
  
  
  plot( x = n_vec, y = unlist( lapply( 1:length(n_vec), FUN = FUN ) ) )
  
  (3 * n_vec + 1) / ( (n_vec + 2) * (n_vec + 3) )
  (3 * n_vec + 1) / ( (n_vec + 2) * (n_vec + 3) ) * (2 * krtsis + n_vec)
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Sampling distributions
  # --------------------------------------------------------------------------
  
  # Exponential population
  beta_param <- 1.1
  x <- rexp( n = 10^7, rate = 1 / beta_param )
  n <- 100
  # z <- scale( x[1:n], TRUE, FALSE ) + beta_param
  z <- x[1:n]
  
  mean(x); beta_param; mean(z)
  var(x); beta_param^2; var(x) * (length(x)-1) / length(x)

  Boot <- boot( data = z,
                statistic = meanBoot,
                R = 10^5 )
  var(Boot$t) * (n-1) / n
  beta_param^2 / n
  var(z) * (n-1) / n / n
  
  
  
  
  # --------------------------------------------------------------------------
  # Expected squared 2-norm FEV-Bias - Bayesian bootstrap
  # --------------------------------------------------------------------------
  
  
  # Finding:
  # Expected squared 2-norm of FEV-biased sample does NOT correspond to expected
  # squared 2-norm of Dirichlet with scaled symmetric alpha': alpha / bias^2
  
  

  betaFUN <- function(idx, alpha) { a <- alpha[idx]; betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  
  bias_vec <- c(1, 3, 5)
  n_vec <- seq(from = 10, to = 10^2, length.out = 10)
  e2n_fev <- matrix( NA, nrow = length(n_vec), ncol = length(bias_vec),
                     dimnames = list(paste0("n=", n_vec), paste0("bias=", bias_vec)) )
  e2n_tilde <- e2n_fev
  e2n_analytic <- e2n_fev
  
  tic <- Sys.time()
  for( i in seq(alon = n_vec) ) {
    
    # Dirichlet parameter
    alpha <- rep(1, n_vec[i])
    wghts <- alpha / n_vec[i]
    
    # Data
    df <- 4
    set.seed(1111)
    x <- rt( n = n_vec[i], df = df )
    z <- x^2
    s_ml <- cov.wt( x = matrix(z, ncol = 1), 
                    wt = wghts, 
                    cor = FALSE, 
                    center = TRUE, 
                    method = "ML" )$cov
    
    for ( j in seq(along = bias_vec) ) {
      
      alpha_tilde <- alpha / bias_vec[j]^2
      beta_moments_tilde <- unlist( lapply( 1:length(alpha_tilde), FUN = betaFUN, 
                                            alpha = alpha_tilde ) )
      e2n_analytic[i, j] <- sum(beta_moments_tilde)
      
      P <- rdirichlet( n = 10^4, alpha = alpha )
      P_tilde <- rdirichlet( n = 10^4, alpha = alpha_tilde )
      e2n_tilde[i, j] <- mean( apply( P_tilde, 1, lpNorm, p = 2 )^2 )
      P_fev <- fevBias( x = P, q = bias_vec[j] )
      e2n_fev[i, j] <- mean( apply( P_fev, 1, lpNorm, p = 2 )^2 )
      
      s_ml * ( 1 + (e2n_analytic[i, j] - 1) / (1 - sum(wghts^2) ) )
      s_ml * ( 1 + (e2n_tilde[i, j] - 1) / (1 - sum(wghts^2) ) )
      s_ml * ( 1 + (e2n_fev[i, j] - 1) / (1 - sum( apply(P_fev, 2, mean)^2 ) ) )
      
      var( apply(P, 1, function(w) { sum(w * z) } ) )
      var( apply(P_tilde, 1, function(w) { sum(w * z) } ) )
      var( apply(P_fev, 1, function(w) { sum(w * z) } ) )
      
    }
  }
  ( toc <- Sys.time() - tic )
  
  1/n_vec
  e2n_analytic
  e2n_tilde
  e2n_fev
  
  
  
  # Does formula for exact standard error still hold?
  
  
  
  
  
  
  
  
  
  
  
  
  plot( as.timeSeries(e2n_analytic), plot.type = "single" )
  
  
  plot( x = n_vec, y = e2n[ ,1]^2 )
  lines( x = n_vec, y = 2 / (n_vec + 1), col = 2 )
  
  plot( x = n_vec, y = e2n_tilde[ ,1]^2 )

  
  
  plot( x = bias_vec, y = e2n[1, ] )
  plot( x = bias_vec, y = e2n[100, ] )
  
  
  colors <- fBasics::divPalette( n = nrow(e2n), "Spectral" )
  plot( x = bias_vec, y = e2n[1, ], ylim = range(e2n) )
  for ( i in 1:nrow(e2n) ) {
    points( x = bias_vec, y = e2n[i, ], col = colors[i] )
  }
  
  colors <- fBasics::divPalette( n = ncol(e2n), "Spectral" )
  plot( x = n_vec, y = e2n[ ,1]^2, ylim = range(e2n^2) )
  for ( i in 1:ncol(e2n) ) {
    points( x = n_vec, y = e2n[ ,i]^2, col = colors[i] )
  }
  lines( x = n_vec, y = 2 / (n_vec + 1), col = 2 )
  
  
  colors <- fBasics::divPalette( n = ncol(e2n_tilde), "Spectral" )
  plot( x = n_vec, y = e2n_tilde[ ,1]^2, ylim = range(na.omit(e2n_tilde)^2) )
  for ( i in 1:ncol(e2n_tilde) ) {
    points( x = n_vec, y = e2n_tilde[ ,i]^2, col = colors[i] )
  }
  lines( x = n_vec, y = 2 / (n_vec + 1), col = 2 )
  lines( x = n_vec, y = e2n_analytic[ ,1], col = 3 )
  lines( x = n_vec, y = e2n_analytic[ ,2], col = 4 )
  
  
  
  
  require(rgl)
  surface3d( x = n_vec/10, y = bias_vec, z = e2n^2 * 100, col = "blue" )
  
  plot3d( x = P[ ,1], y = P[ ,2], z = P[ ,3] )
  
  
  
  
  
    