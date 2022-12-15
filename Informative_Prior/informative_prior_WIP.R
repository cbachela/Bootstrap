  
  
  ############################################################################
  ### INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     04.06.2021
  # First version:    04.06.2021
  # --------------------------------------------------------------------------
  
  
  # Research question:
  # Can we find an exact formulation for the standard error when sampling 
  # probabilities are not all equal ?
  
  
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
  # Toy example
  # --------------------------------------------------------------------------
  
  n <- 6
  B <- 10^4
  n_boot <- 50
  
  df <- 5
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  
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
  ent_parts <- apply( parts_mat / n, 2, entropy, exponential = FALSE, base = exp(1) )
  ent_parts <- apply( parts_mat / n, 2, entropy, exponential = TRUE, base = exp(1) )
  
  
  

  parts_mat
  draws
  
  
  
  param_multinom <- (1:n) / sum(1:n)
  # param_multinom <- rep(1/n, n)
  param_multinom
  prob_prior <- apply( parts_mat, 2, dmultinom, prob = param_multinom, log = FALSE )
  prob_prior
  apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  prob
  H <- -sum( param_multinom * log(param_multinom, base = 2) )
  
  
  
  # Need to enumerate all permutations within each level set

  # # ?gtools::permutations
  # # gtools::permutations( n = n, r = n, v = c(2, 1, 0), repeats.allowed = FALSE )
  # # gtools::permutations( n = n, r = n, v = c(2, 1, 0), repeats.allowed = TRUE )
  # 
  # lP <- list()
  # for( j in 1:ncol(parts_mat) ) {
  #   perm <- try( gtools::permutations( n = n, r = n, v = parts_mat[ ,j], 
  #                                      repeats.allowed = FALSE ) )
  #   if ( inherits(perm, "try-error") ) {
  #     lP[[j]] <- matrix( parts_mat[ ,j], nrow = 1 )
  #   } else {
  #     lP[[j]] <- perm
  #   }
  # }
  # lP
  # 
  # lpP <- lapply( lP, FUN = function(p) { apply( p / n, 1, lpNorm, p = 2) } )
  # lpP
  
  
  FUN <- function(i)
  {
    P <- RP::permutations(parts_mat[ ,i])
    P <- matrix( P[!duplicated(P), ], ncol = n )
  }
  tic <- Sys.time()
  lP <- lapply( 1:ncol(parts_mat), FUN = FUN )
  (toc <- Sys.time() - tic)
  # lP
  
  # # Alternatively:
  # 
  # N <- 10^5
  # P <- matrix( NA, nrow = N, ncol = nrow(parts_mat) )
  # lP <- list()
  # lTheta <- list()
  # for ( k in 1:ncol(parts_mat) ) {
  #   for ( j in 1:N ) {
  #     P[j, ] <- sample( x = parts_mat[ ,k], replace = FALSE )
  #   }  
  #   lP[[k]] <- matrix( P[!duplicated(P), ], ncol = n )
  #   lTheta[[k]] <- (lP[[k]] / n) %*% z
  # }
  # cbind( draws[2, ], lapply( lP, nrow ) )
  
  
  # Theta
  lTheta <- lapply( lP, FUN = function(P) { (P / n) %*% z } )
  
  # Variance of lTheta's
  theta_var <- unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(theta))^2 ) / length(theta) } ) )

  (lp_parts^2 - 1/n) * var(z)
  theta_var
  
  
  
  
  
  # Probability weighted variance of lTheta's
  
  
  lprob <- lapply( lP, function(P) { apply( P, 1, dmultinom, prob = param_multinom, log = FALSE ) } )
  sum( unlist( lprob ) )
  prob_tot_w <- unlist( lapply( lprob, sum ) )
  
  
  # Covariance under exponentially weighted probability measure
  cov.prob <- function(X, p)
  {
    mu <- t(X) %*% p
    scnd_mom <- t(X) %*% (X * (p %*% matrix(1, 1, ncol(X))) )
    scnd_mom <- (scnd_mom + t(scnd_mom)) / 2
    sigma <- scnd_mom - mu %*% t(mu)
    return( sigma )
  }
  
  weightedVariance <- function(x, p, unbiased = FALSE) 
  { 
    ans <- sum( p * (x - sum(p * x))^2 ) 
    if ( isTRUE(unbiased) ) {
      ans <- ans * n / (n - 1)
    }
    return( ans )
  }
  
  cov.prob( X = matrix(z, ncol = 1), p = param_multinom )
  weightedVariance( x = z, p = param_multinom, unbiased = FALSE )
  cov.wt( x = matrix(z, ncol = 1), wt = param_multinom, cor = FALSE, center = TRUE, method = "ML" )
  sum( z * param_multinom )
  
  var(z)
  cov.wt( x = matrix(z, ncol = 1), wt = param_multinom, cor = FALSE, center = TRUE, method = "unbiased" )
  weightedVariance( x = z, p = param_multinom, unbiased = TRUE )
 
  
  # idx <- length(lP)
  idx <- 1
  P <- lP[[idx]]
  prob_l <- lprob[[idx]]
  prob_l_scl <- prob_l / sum(prob_l)
  lp2 <- apply( P / n, 1, lpNorm, p = 2 )
  
  
  cov.prob( X = lTheta[[idx]], p = prob_l_scl )
  weightedVariance( x = lTheta[[idx]], p = prob_l_scl )
  
  cov.prob( X = matrix(z, ncol = 1), p = param_multinom )
  weightedVariance( x = z, p = param_multinom )
  
  
  
  
  
  
  ( sum(lp2^2 * prob_l_scl) - min(lp2^2 * prob_l_scl) ) * var(z)
  ( sum(lp2^2 * prob_l_scl) - lpNorm(param_multinom, p = 2)^2 ) * var(z)
  ( sum(lp2^2 * prob_l_scl) - lpNorm(param_multinom, p = 2)^2 ) * var(z)
  
  
  ( sum(lp2^2 * prob_l_scl) - 1/n ) * cov.prob( matrix(z, ncol = 1), param_multinom )
  ( sum(lp2^2 * prob_l_scl) - min(lp2^2 * prob_l_scl) ) * cov.prob( matrix(z, ncol = 1), param_multinom )
  ( sum(lp2^2 * prob_l_scl) - min(lp2^2 * prob_l_scl) ) * weightedVariance( z, param_multinom )
  
  ( sum(lp2^2 * prob_l_scl) - lpNorm(param_multinom, p = 2)^2 ) * weightedVariance( z, param_multinom )
  
  
  rhs <- cov.prob( X = lTheta[[idx]], p = prob_l_scl )
  rhs
  

  # Unweighted and weighted mean
  mean(z); mean(unlist(lTheta)) # same same
  sum(z * param_multinom); sum( unlist(lTheta) * unlist(lprob) ) # same same
  a <- unlist( lapply( 1:length(lTheta), FUN = function(i) { sum(lTheta[[i]] * lprob[[i]]) } ) )
  sum(a) # same same
  b <- unlist( lapply( 1:length(lTheta), FUN = function(i) { sum(lTheta[[i]] * lprob[[i]] / sum(lprob[[i]])) } ) )
  sum( b * unlist( lapply( lprob, sum ) ) )  # same same
  
  
  
  # Weighted variance
  v <- unlist( lapply( 1:length(lTheta), FUN = function(i) {  sum( lprob[[i]] * (lTheta[[i]] - sum(a))^2 ) } ) )
  v
  weightedVariance( x = unlist(lTheta), p = unlist(lprob), unbiased = FALSE ) # same same
  sum(v) # same same
  v2 <- unlist( lapply( 1:length(lTheta), FUN = function(i) {  sum( lprob[[i]] / sum(lprob[[i]]) * (lTheta[[i]] - sum(a))^2 ) } ) )
  v2
  sum( v2 * unlist( lapply( lprob, sum ) ) )  # same same
  v3 <- unlist( lapply( 1:length(lTheta), FUN = function(i) {  sum( lprob[[i]] / sum(lprob[[i]]) * (lTheta[[i]] - b[i])^2 ) } ) )
  sum( v3 * unlist( lapply( lprob, sum ) ) )
  
  

  plot( x = lp_parts^2, y = v2 )
  plot( x = lp_parts^2, y = v3 )
  plot( x = exp(ent_parts), y = v2 )
  plot( x = lp_parts^2, y = unlist(lapply(lTheta, var)), type = "o" ) # ??
  plot( x = lp_parts^2, y = unlist(lapply(lTheta, function(theta) { sum( (theta - sum(a))^2) / length(theta) })), type = "o" )
  
  
  
  # Typicality
  
  H <- entropy( p = param_multinom, exponential = FALSE, base = 2 )
  H <- -sum( param_multinom * log(param_multinom, base = 2) )  # same same
  eps <- 1e-03
  lb <- 2^(n * (-H + eps))
  ub <- 2^(n * (-H - eps))
  lb; ub
  lPent <- lapply( lP, function(P) { apply( P, 1, entropy, exponential = FALSE, base = 2 ) } )  
  # lPent
  
  # AEP
  range( -1/n * log(unlist(lprob), base = 2) )   # --> to H(X)
  H
  boxplot( -1/n * log(unlist(lprob), base = 2) )
  abline( h = H )
  
  
  1 / H * weightedVariance(x = z, p = param_multinom)
  
  
  
    
  
  theta_mu <- unlist( lapply( 1:length(lTheta), FUN = function(i) { mean(lTheta[[i]]) } ) )
  theta_mu_w <- unlist( lapply( 1:length(lTheta), FUN = function(i) { sum( lprob[[i]] * lTheta[[i]] ) } ) )
  barplot( t(cbind(theta_mu, theta_mu_w)), beside = TRUE, col = 1:2 )
  
  theta_var_w <- unlist( lapply( 1:length(lTheta), FUN = function(i) {  cov.prob( X = matrix(lTheta[[i]], ncol = 1), p = lprob[[i]] / sum(lprob[[i]]) ) } ) )
  barplot( t(cbind(theta_var, theta_var_w)), beside = TRUE, col = 1:2 )
  
  
  
  x <- lp_parts^2
  cbind( x, theta_var_w )
  plot( x = x, y = theta_var_w )
  plot( x = theta_var, y = theta_var_w )
  plot( x = ent_parts, y = theta_var_w )
  
  
  
  # Expected 2-norm under prior probabilities
  lp_lP <- unlist( lapply( lP, function(P) { apply( P / n, 1, lpNorm, p = 2 ) } ) )
  elp2 <- sum( lp_lP^2 * unlist(lprob) )
  elp2
  
  
  # Ordinary bootstrap
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_w <- boot_mat
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        weights = param_multinom,
                        R = B )
    boot_mat_w[ ,j] <- Boot$t
    
  }
  mean( apply( boot_mat, 2, var ) )
  ( sum( lp_parts^2 * prob_tot ) - 1/n ) * var(z)
  
  Vz <- cov.wt( x = matrix(z, ncol = 1), wt = param_multinom, cor = FALSE, center = TRUE, method = "unbiased" )
  ( elp2 - lpNorm(param_multinom, p = 2)^2 ) * Vz$cov
  mean( apply( boot_mat_w, 2, var ) )
  
  
  
  
  
  


  
  
  
  
  
  ################
  
  
 
  draws2prob
  a <- draws[1, ]
  b <- draws[2, ]
  d <- as.numeric( a %*% b )
  p <- a / d  
  p  
  
  prob_prior
  
  (1/n)^n
  
  
  d2 <- a / prob_prior
  p2 <- a / d2
  p2
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact BB standard error based on expected 2-norm - with informative priors
  # --------------------------------------------------------------------------
  
  n <- 80
  n_boot <- 50
  B <- 10^5
  alpha <- (1:n) / sum(1:n) * n
  # alpha <- rep(1, n)
  alpha_std <- alpha / sum(alpha)
  
  df <- 5
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  # Bayesian bootstrap
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP_Dir <- list()
  lNorm <- list()
  for ( j in 1:n_boot ) {
    lP_Dir[[j]] <- rdirichlet( n = B, alpha = alpha )
    lNorm[[j]] <- apply( lP_Dir[[j]], 1, lpNorm, p = 2 )
    bb_mat[ ,j] <- apply( lP_Dir[[j]], 1, function(p) { sum(p * z) } )
  }

  
  # Expected 2-norm
  expected_2normSq <- unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) )
  
  
  # Analytic solution
  FUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = FUN ) )
  beta_moments  
  sum(beta_moments)
  2 / (n + 1)
  
  expected_2normSq
  sum(beta_moments)
  mean(expected_2normSq)
  
  Vz <- cov.wt( x = matrix(z, ncol = 1), wt = alpha_std, cor = FALSE, center = TRUE, method = "unbiased" )
  stderr_exact_bb <- sqrt( Vz$cov * ( sum(beta_moments) - lpNorm(alpha_std, p = 2)^2 ) )
  stderr_bb <- apply( bb_mat, 2, sd )

  stderr_exact_bb  
  stderr_bb       # different
  mean(stderr_bb)
  
  boxplot( stderr_bb )
  abline( h = stderr_exact_bb )
  
 
 
  
  
  
  # --------------------------------------------------------------------------
  # Typical set
  # --------------------------------------------------------------------------
  
  H <- -sum( param_multinom * log(param_multinom, base = 2) )
  
  
  eps <- 1e-03
  a <- 2^(n * (-H + eps))
  b <- 2^(n * (-H - eps))
  boxplot( log(prob_tot_w) )
  abline( h = log(a) )
  abline( h = log(b), col = 2 )
    
  range( -1/n * log(prob_tot_w, base = 2) )
  H  
  
  
  
  
  
  
  
  
  
  # # --------------------------------------------------------------------------
  # # Permutations
  # # --------------------------------------------------------------------------
  # 
  # permutations <- function( x ) 
  # {
  #   if ( length(x) == 1 ) {
  #     return( x )
  #   }
  #   else {
  #     tmp <- table(x)
  #     if ( length(x) == sum(x) ) {
  #       num <- factorial( sum(x) )
  #     } else {
  #       num <- factorial( length(x) )
  #     }
  #     n_row <- num / prod( factorial( tmp ) )
  #     # res <- matrix( x, nrow = n_row, ncol = length(x), byrow = TRUE )
  #     res <- matrix( 0, nrow = n_row, ncol = length(x) )
  #     res[1, ] <- x
  #     i <- 2
  #     for ( k1 in seq(along = x) ) {
  #       for ( k2 in seq(along = x) ) {
  #         if ( k1 != k2 ) {
  #           y <- x
  #           y[k1] <- x[k2]
  #           y[k2] <- x[k1]
  #           if ( !any( apply( res[1:i, ], 1, function(x) { all(x == y) } ) ) ) {
  #             res[i, ] <- y
  #             i <- i + 1
  #           }
  #         }
  #       }
  #     }
  #   }
  #   return( res )
  # }
  # 
  # 
  # 
  # debugonce( permutations )
  # permutations( x = c(2, 1, 0) )
  
  
  
  