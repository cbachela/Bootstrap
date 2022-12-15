  
  
  ############################################################################
  ### EXACT SAMPLE DISTRIBUTION OF THE MEAN
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     14.04.2022
  # First version:    14.04.2022
  # --------------------------------------------------------------------------

  
  
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  # Check if a^{\top} B a = \sum_i a_i^2 b_i^2 + \sum_i \sum_{j, i \noteq j} a_i a_j b_i b_j
  # --> it does!
  n <- 10
  set.seed(1111)
  a <- runif(n)
  b <- runif(n)
  
  t(a) %*% (b %*% t(b)) %*% a
  sum( a^2 * b^2 ) + sum( (b %*% t(b))[upper.tri(b %*% t(b))] * (a %*% t(a))[upper.tri(a %*% t(a))] ) * 2
  
  
  mu <- 0.1
  a <- rnorm(10^7, mu, 1)
  b <- rchisq(10^7, df = 1)
  y <- a^2
  mean(y)
  
  
  
  # Sections:
  # Bayesian bootstrap
  # Moments of the approximate boostrap distribution based on norms
  # When is \sum_{i=1}^n (x_i - \bar{x}_i)^2 = ||x||_2^2 - ||\bar{x}||_2^2 ?
  # Scalefactor surface
  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
  n <- 3
  
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
  
  
  
  # --------------------------------------------------------------------------
  # Permutations
  # --------------------------------------------------------------------------
  
  lP <- list()
  lTheta <- list()
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
  gridnodes_all <- do.call( rbind, lP )
  Grid <- gridnodes_all / n
  Grid
  dim(Grid)
  comb( n = 2*n-1, k = n )
  
  lProb <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  prob_k <- unlist(lProb)
 
  
  ###
  # prob_k <- 1:nrow(Grid) / sum(1:nrow(Grid))
  ###
  
  
  
  # --------------------------------------------------------------------------
  # Normal data
  # --------------------------------------------------------------------------
  
  m1 <- 0.1 # population mean
  m2 <- 1   # population variance
  set.seed(1:10)
  x <- rnorm( 10^7, m1, sqrt(m2) )
  
  # Take a sample
  z <- x[1:n]
  
  # Sample mean
  z_bar <- mean(z)
  m2_hat <- var(z)
  m2_hat_ml <- m2_hat * (n-1) / n
    
  # Theoretical variance of z_bar
  m2 / n          # using population variance
  m2_hat / n      # using sample variace
  m2_hat_ml / n   # using biased sample variance
  
  
  # Bootstrap (classical)
  n_boot <- 50
  B <- 10^4 + 1
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  tic <- Sys.time()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  apply(boot_mat, 2, var)

  
  
  
  m2 / n
  m2_hat / n
  range( apply( boot_mat, 2, var ) )
  m2_hat_ml / n; mean( apply( boot_mat, 2, var ) )  # close enough
  
  

  
  
  # Using (grid) weights
  
  ewi2 <- 2/n - 1/n^2
  sum( apply( Grid, 1, function(g) { sum( g^2 ) } ) * prob_k )
  sum( prob_tot * lp_parts^2 ) # same same
  ewi2  # same same
  
  Grid

  # Expectation under the distribution over the Grid
  # \mathds{E}_G (\theta - \mathds{E}_G(\theta))^2
  b <- z - mean(z)
  Ans <- matrix(0, nrow = nrow(Grid), ncol = 2)
  for ( k in 1:nrow(Grid) ) {
    g <- as.numeric(Grid[k, ])
    Ans[k, 1] <- sum(g^2 * b^2)
    tmp <- 0
    for ( i in 1:(n-1) ) {
      for ( j in (i+1):n ) {
        tmp <- tmp + g[i] * g[j] * b[i] * b[j] * 2
      }
    }
    Ans[k, 2] <- tmp
  }
  ans <- apply( Ans, 1, sum )
  
  apply( Ans, 2, function(a) { sum(a * prob_k) } )
  sum( apply( Ans, 2, function(a) { sum(a * prob_k) } ) )
  sum( prob_k * ans ); m2_hat_ml / n # same same
  Bmat <- b %*% t(b)
  ans2 <- apply( Grid, 1, function(g) { t(g) %*% Bmat %*% g } )
  sum( ans2 * prob_k ) # same same
  range( ans - ans2 )
  
  # Alternatively, use expected weights and expected squared weights
  wsq_bar <- apply( Grid^2 * prob_k, 2, sum )
  sum( wsq_bar * b^2 ); sum( Ans[ ,1] * prob_k )   # same same
  
  
  w_bar <- apply( Grid * prob_k, 2, sum )
  t(w_bar) %*% Bmat %*% w_bar - sum( w_bar^2 * diag(Bmat) )
  sum( Ans[ ,2] * prob_k ) # different
  
  tmp <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp <- tmp + w_bar[i] * w_bar[j] * b[i] * b[j] * 2
    }
  }
  tmp
  
  
  # Or use expected product w_i w_j
  W <- t( apply( Grid, 1, function(g) { (g %*% t(g))[upper.tri(g %*% t(b))] } ) )
  wiwj_bar <- apply( W * prob_k, 2, sum )
  wisq_bar <- apply( Grid^2, 2, mean )
  wiwj_bar
  wsq_bar
  
  wiwj_bar_mat <- matrix(NA, nrow = n, ncol = n, byrow = TRUE)
  wiwj_bar_mat[upper.tri(wiwj_bar_mat)] <- wiwj_bar
  mat <- t(wiwj_bar_mat)
  mat[upper.tri(mat)] <- wiwj_bar
  wiwj_bar_mat <- t(mat)
  wiwj_bar_mat
  tmp <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp <- tmp + wiwj_bar_mat[i, j] * b[i] * b[j] * 2
    }
  }
  tmp; sum( Ans[ ,2] * prob_k )   # same same
  sum(wsq_bar * b^2) + tmp; m2_hat_ml / n  # same same
  
  

  # Also same same:
  w <- w_bar <- rep(1/n, n)
  sum( w_bar * (z - z_bar)^2 ) * 1 / (1 - sum(w_bar^2)) * ( ewi2 - sum(w_bar^2) )
  m2_hat_ml / n
  
  
  # Only consider partitions
  ans3 <- apply( t(parts_mat) / n, 1, function(g) { t(g) %*% Bmat %*% g } )
  sum( ans3 * prob_tot) # different
  
  

  
  # Compare to case where weights are fixed, data are random z \sim \hat{F}
  sum(w_bar^2) * m2_hat_ml
  
  
  
  
  
  
  
  
  
  ######
  
  # Assuming weights are fixed, data are random
  w <- rep(1/n, n)
  w_bar <- w
  m2 * sum(w^2);  m2 / n
  m2_hat * sum(w^2);  m2_hat / n
  
  b <- z - z_bar
  
  
  vec1 <- NULL
  vec2 <- NULL
  for( i in 1:n ) {
    # vec1 <- c(vec1, w_bar[i]^2 * b[i]^2 )
    vec1 <- c(vec1, w_bar[i]^2 * m2_hat )
    for( j in 1:n ) {
      if ( i != j ) {
        vec2 <- c(vec2, w_bar[i] * w_bar[j] * b[i] * b[j])
      } 
    }
  }
  vec1
  vec2
  sum(vec1) + sum(vec2)
  
  sum(vec1); sum( w^2 * m2_hat ) # same same
  range( apply( boot_mat, 2, var ) ); mean( apply( boot_mat, 2, var ) )
  
  
  
  ###
  B <- b %*% t(b)
  t(w_bar) %*% B %*% w_bar  #// variance is zero bcs. weights are fixed and data are fixed
  ###
  
  
  
  
  
  # Assuming data are fixed, weights are random
  
  w <- rep(1/n, n)
  w_bar <- rep(1/n, n)
  # ewi2 <- 2/n^2 - 1/n^3
  ewi2 <- 2/n - 1/n^2
  
  sum(prob_tot * lp_parts^2); ewi2  # same same
  
  
  vec1 <- NULL
  vec2 <- NULL
  for( i in 1:n ) {
    vec1 <- c(vec1, ewi2 / n * b[i]^2 )
    for( j in 1:n ) {
      if ( i != j ) {
        vec2 <- c(vec2, w_bar[i] * w_bar[j] * b[i] * b[j])
      } 
    }
  }
  vec1; ewi2 / n * b^2
  sort(vec2); sort( rep( (b %*% t(b))[upper.tri(b %*% t(b))] * 1/n^2, 2 ) )
  sum(vec1) + sum(vec2)
  sum( w^2 * m2_hat_ml )
  sum(vec1)
  sum(vec2)
  
  m2 * 1/n; m2 * sum(w^2)
  m2_hat / n; m2_hat * sum(w^2)
  m2_hat_ml / n; m2_hat_ml * sum(w^2)
  
  range( apply( boot_mat, 2, var ) ); mean( apply( boot_mat, 2, var) )
  sum(vec1) + sum(vec2)

  
  # Exact ordinary bootstrap
  Vz_unb <- var(z)
  a <- sum( (z - mean(z))^2 )
  Vz_ml <- a / n
  stderr_parts <- sqrt( var(z) * ( sum(prob_tot * lp_parts^2) - 1/n ) )
  stderr_analytic <- sqrt( Vz_ml / n )
  
  stderr_analytic^2; m2_hat_ml / n; stderr_parts^2  # same same same
  
 
  
  
  
  #####################
  
  
  
  
  
  n <- 3
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  w_bar <- rep(1/n, n)
  m2 <- sum( w_bar * (z - mean(z))^2 )  
  
  
  m2 * 1 / (1 - sum(w_bar^2)) * ( 2/n - 1/n^2 - sum(w_bar^2) )
  m2 * sum(w_bar^2)
  m2 * 1/n  
  2 * m2 - m2 / n
  
  
  
  b <- z
  w_bar[1]^2 * z[1]^2 + w_bar[2]^2 * z[2]^2 + 2 * w_bar[1] * w_bar[2] * z[1] * z[2] - sum( mean(z) * w_bar * z )
  
  
  b <- z - mean(z)
  A <- b %*% t(b)
  A
  
  t(w_bar) %*% A %*% w_bar
  sum( w_bar * (z - mean(z)) )^2
  w_bar[1]^2 * b[1]^2 + w_bar[2]^2 * b[2]^2 + 2 * w_bar[1] * w_bar[2] * b[1] * b[2]
  
  
  ewi2 <- 2/n^2 - 1/n^3
  b <- z - mean(z)
  ewi2 * b[1]^2 + ewi2 * b[2]^2 + 2 * w_bar[1] * w_bar[2] * b[1] * b[2]
  (ewi2 * b[1]^2 + ewi2 * b[2]^2 + 2 * w_bar[1] * w_bar[2] * b[1] * b[2]) * 2
  
  
  ( sum( ewi2 * b^2 ) + 
      2 * prod( b[c(1, 2)] ) * prod( w_bar[c(1, 2)] ) + 
      2 * prod( b[c(1, 3)] ) * prod( w_bar[c(1, 3)] ) +
      2 * prod( b[c(3, 2)] ) * prod( w_bar[c(3, 2)] ) )
  
  m2 * 1/n  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # --------------------------------------------------------------------------
  
  # Generate Normally distributed data
  m1 <- 0.1 # population mean
  m2 <- 1   # population variance
  set.seed(1:10)
  x <- rnorm( 10^7, m1, sqrt(m2) )
  
  # Take a sample
  z <- x[1:n]
  z_bar <- mean(z)
  m2_hat_ml <- sum( (z - mean(z))^2 ) / n
  
  
  # Run sampler
  w_bar <- rep(1/n, n)
  w <- w_bar
  alpha <- w_bar * n
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = alpha )
    lP[[j]][is.na(lP[[j]])] <- 0
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  
  # Empirical expected squared 2-norm
  l_ewi2_hat <- lapply( lP, FUN = function(P) { apply(P, 1, function(p) { sum(p^2) } ) } )
  ewi2_hat <- mean( unlist( l_ewi2_hat ) )
  
  # Empirical expected E(w_i w_j)
  l_ewiwj_hat <- lapply( lP, FUN = function(P) { apply( P, 1, function(p) { pmat <- p %*% t(p); pmat[upper.tri(pmat)] } ) } )
  ewiwj_mat_hat <- do.call( cbind, lapply( l_ewiwj_hat, FUN = function(x) { apply(x, 1, mean) } ) )
  apply( ewiwj_mat_hat, 1, mean )
  ewiwj_hat <- mean( ewiwj_mat_hat )
  
  
  # Theoretical E(w_i^2)
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  beta_moments; 2 / (n * (n + 1)); ewi2_hat / n
  ew1sq <- beta_moments[1]
  
  # Theoretical E(w_i w_j)
  ewiwj <- (alpha[i] * alpha[j]) / (sum(alpha) * (sum(alpha) + 1))
  ewiwj
  weiwj_hat
  
  # Compute second moment 
  b <- (z - z_bar)
  tmp <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp <- tmp + b[i] * b[j]
    }
  }
  
  ewi2_hat/n * sum( (z - z_bar)^2 ) + 2 * ewiwj_hat * tmp
  ew1sq * sum( (z - z_bar)^2 ) + 2 * ewiwj * tmp
  M2Analytic( z = z, method = "Bayesian" ) # same same
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Moments of the approximate boostrap distribution based on norms
  # --------------------------------------------------------------------------
  
  # Question: 
  # Can the moments of the approximate bootstrap distributions be computed by 
  # m_2 ( \frac{1}{B} \sum_{i=1}^B (\omega_i^2)^b - ||\frac{1}{B} \sum_{i=1}^B (\omega_2)^b)^2||_2^2 ) ?
  
  
  # Parameters
  n <- 3
  B <- 21
  
  # Partitions and permutations
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
  
  # Compute permutations
  lP <- list()
  lTheta <- list()
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
  gridnodes_all <- do.call( rbind, lP )
  Grid <- gridnodes_all / n
  Grid
  dim(Grid)
  comb( n = 2*n-1, k = n )
  
  lProb <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  prob_k <- unlist(lProb)
  
  
  # Generate Normally distributed data
  m1 <- 0.1 # population mean
  m2 <- 1   # population variance
  set.seed(1:10)
  x <- rnorm( 10^7, m1, sqrt(m2) )
  
  # Take a sample
  z <- x[1:n]
  z_bar <- mean(z)
  m2_hat_ml <- sum( (z - mean(z))^2 ) / n
  
  
  # Classical bootstrap by hand
  idx <- numeric(B)
  for ( i in 1:B ) { 
    idx[i] <- sample( 1:nrow(Grid), size = 1, replace = TRUE, prob = prob_k )
  }
  theta <- apply( Grid[idx, ], 1, function(g) { sum(g * z) } )
  
  # 2nd moment of the bootstrap theta's
  mu_2_B <- sum( (theta - mean(theta))^2 ) / length(theta)
  
  # Expected (2-norm)^2 of the sampled weights
  exp_norm_2_sq <- mean( apply( Grid[idx, ], 1, function(g) { sum(g^2) } ) )
  exp_norm_2_sq; 2/n - 1/n^2 # close
  
  # Average weights vector
  w_bar_emp <- apply( Grid[idx, ], 2, mean )
  # Weighted average z
  z_bar_emp <- sum(w_bar_emp * z)
  mean(theta); z_bar; z_bar
  
  # Exact bootstrap 2nd moment
  mu_2_exact <- M2Exact( z = z, w_bar = w_bar, exp2norm2 = 2/n - 1/n^2 )
  mu_2_exact
  m2_hat_ml / n
  mu_2_B
  
  
  
  b <- z - z_bar_emp # notice the use of _emp
  Ans <- matrix(0, nrow = length(idx), ncol = 2)
  for ( k in 1:length(idx) ) {
    g <- as.numeric(Grid[idx[k], ])
    Ans[k, 1] <- sum(g^2 * b^2)
    tmp <- 0
    for ( i in 1:(n-1) ) {
      for ( j in (i+1):n ) {
        tmp <- tmp + g[i] * g[j] * b[i] * b[j] * 2
      }
    }
    Ans[k, 2] <- tmp
  }
  ans <- apply( Ans, 1, sum )
  
  apply( Ans, 2, mean )
  sum( apply( Ans, 2, mean) ); mean(ans); mu_2_B  # same same
 
  
  
  # Can we use expected weights? --> Yes!
  
  wsq_bar <- apply( Grid[idx, ]^2, 2, mean )
  sum( wsq_bar * b^2 ); mean( Ans[ ,1] )   # same same
  
  # Expected product w_i w_j
  W <- t( apply( Grid[idx, ], 1, function(g) { (g %*% t(g))[upper.tri(g %*% t(b))] } ) )
  wiwj_bar <- apply( W, 2, mean )
  wisq_bar <- apply( Grid[idx, ]^2, 2, mean )
  wiwj_bar_mat <- matrix(NA, nrow = n, ncol = n, byrow = TRUE)
  wiwj_bar_mat[upper.tri(wiwj_bar_mat)] <- wiwj_bar
  mat <- t(wiwj_bar_mat)
  mat[upper.tri(mat)] <- wiwj_bar
  wiwj_bar_mat <- t(mat)
  wiwj_bar_mat
  tmp <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp <- tmp + wiwj_bar_mat[i, j] * b[i] * b[j] * 2
    }
  }
  tmp; mean( Ans[ ,2] )   # same same
  sum(wsq_bar * b^2) + tmp; mean(ans); mu_2_B  # same same
  
  
  
  # debugonce( M2Exact )
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = exp_norm_2_sq )
  #
  M2Exact( z = z, w_bar = w_bar_emp, exp2norm2 = exp_norm_2_sq )
  M2Exact( m2 = sum( w_bar_emp * (z - z_bar_emp)^2 ) / (1 - sum(w_bar_emp^2)),
           w_bar = w_bar_emp, 
           exp2norm2 = exp_norm_2_sq )
  #
  M2Exact( m2 = sum( w_bar_emp * (z - z_bar_emp)^2 ) / (1 - sum(w_bar_emp^2)),
           w_bar = w_bar, 
           exp2norm2 = 2/n - 1/n^2 )

  
  
  
  
  
  # --------------------------------------------------------------------------
  # When is \sum_{i=1}^n (x_i - \bar{x}_i)^2 = ||x||_2^2 - ||\bar{x}||_2^2 ?
  # --------------------------------------------------------------------------
  
  # Answer: holds for all elements x on the simplex.
  
  
  n_sim <- 10^4
  set.seed(1234)
  n_vec <- seq(from = 2, to = 10^3, length.out = 30)
  lWbar <- list()
  delta_vec <- n_vec * NA
  
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    alpha <- runif(n = n)
    alpha <- alpha / sum(alpha)
    P <- rdirichlet( n = n_sim, alpha = alpha )
    w_bar <- apply(P, 2, mean)
    lpsq_w_bar <- sum(w_bar^2)
    
    a <- apply( P, 1, function(p) { sum( (p - w_bar)^2 ) } )
    b <- apply( P, 1, function(p) { sum( p^2 ) - lpsq_w_bar } )
    
    tmp <- cbind( apply( P, 1, function(p) { sum( p * w_bar ) } ),
                  apply( P, 1, function(p) { sum( w_bar^2 ) } ) )
    apply( tmp, 2, mean )
    
    w_bar[1]
    1/n
    sum(w_bar^2)
    
    
    delta_vec[i] <- mean(a) - mean(b)
    lWbar[[i]] <- w_bar
    
  }
  
  plot( delta_vec )
  
  
  
  
  
  