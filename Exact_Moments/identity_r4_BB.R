  
  
  ############################################################################
  ### IDENTITY FOR r=3, DIRICHLET CASE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.07.2022
  # First version:    20.07.2022
  # --------------------------------------------------------------------------
  
  
  # Description:
  # This scripts defines an (arbitrarily) real-valued vector z and computes the 
  # RHS and LHS of the identity for the case:
  # - Weights follow a Dirichlet distribution. I
  
  
  
  
  
  require(partitions)
  require(gtools)
  
  
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Data and parameters
  # --------------------------------------------------------------------------
  
  # Define real-valued vector z
  n <- 4
  set.seed(1111)
  z <- rt(n = n, df = 3)^2
  
  w_bar <- rep(1/n, n)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar
  
  # Number of simulations (for the Dirichlet sampling case)
  n_sim <- 10^6
  
  
 
  

  
  # --------------------------------------------------------------------------
  # Weights follow a flat Dirichlet distribution
  # --------------------------------------------------------------------------
  
  # Sample weights from a Dirichlet distribution
  alpha <- rep(1, n) 
  w_bar <- alpha / sum(alpha)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  
  # Compute theta_hat for all sampled weights
  theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
  
  # Norms on weights
  expected_2norm2 <- mean( apply( samples, 1, function(w) { sum(w^2) } ) )
  expected_3norm3 <- mean( apply( samples, 1, function(w) { sum(w^3) } ) )
  exact_expected_2norm2 <- (1 + lambda) / (1 + lambda * n)
  # exact_expected_3norm3 <- ((1 + lambda) * (2 + lambda)) / ((1 + n * lambda) * (2 + n * lambda))
  exact_expected_3norm3 <- (lambda^2 + 3 * lambda + 2) / ( (n * lambda)^2 + 3 * n * lambda + 2 )
  exact_expected_3norm3
  
  expected_2norm2; exact_expected_2norm2  # should be similar
  expected_3norm3; exact_expected_3norm3  # should be similar
  
  
  # Compute RHS of identity
  RHS <- M3Exact( z = z,
                  w_bar = w_bar,
                  exp3norm3 = exact_expected_3norm3,
                  # exp2norm2 = exact_expected_2norm2,
                  expw2 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  RHS  
  
  
  
  # Compute LHS using E(w_i^3), E(w_i^2 w_j) and E(w_i w_j, w_k)
  alpha0 <- sum(alpha)
  
  wi3_bar <- unlist( lapply( alpha, FUN = betaMoments, a0 = sum(alpha), m = 3 ) )
  ( (alpha[1] + 2) * (alpha[1] + 1) * alpha[1] ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  wi3_bar
  
  w1w2 <- mean( apply( samples, 1, function(w) { prod(w[1:2] ) } ) )
  ( alpha[1] * alpha[2] ) / (alpha0 * (alpha0 + 1) )
  w1w2 # same same
  
 
  
  
  # E(w_i w_j, w_k)
  
  # Empirical
  FUN <- function(g)
  {
    c( prod(g[c(1, 2, 3)]), 
       prod(g[c(1, 2, 4)]),
       prod(g[c(1, 3, 4)]),
       prod(g[c(2, 3, 4)]) )
  }
  W <- t( apply( samples, 1, FUN ) )
  wiwjwk_bar<- apply( W, 2, mean )
 
  # Exact
  exactEwiwjwk <- function(alpha)
  {
    alpha0 <- sum(alpha)
    A <- cbind( alpha[c(1, 2, 3)],
                alpha[c(1, 2, 4)],
                alpha[c(1, 3, 4)],
                alpha[c(2, 3, 4)] )
    denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- apply( A, 2, function(a) { prod(a) / denom } )
    return( ans )
  }
  
  exactEwiwjwk( alpha = alpha )
  wiwjwk_bar # same same
  
  
  
  # E(w_i^2 w_j)
  
  # Empirical
  W <- t( apply( samples, 1, function(g) { (g^2 %*% t(g))[upper.tri(g %*% t(g))] } ) )
  if ( nrow(W) == 1 ) W <- t(W)
  wi2wj_bar <- apply( W, 2, mean )
 
  # Exact
  exactEwi2wj <- function( alpha ) 
  {
    alpha0 <- sum(alpha)
    A <- cbind( alpha[c(1, 2)],
                alpha[c(1, 3)],
                alpha[c(1, 4)],
                alpha[c(2, 3)],
                alpha[c(2, 4)],
                alpha[c(3, 4)] )
    denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- apply( A, 2, function(a) { ( (a[1] + 1) * a[2] ) / denom } )
    return( ans )
  }
  
  exactEwi2wj( alpha = alpha )
  wi2wj_bar
  
  
  
  
  
  
  
  

  
  # --------------------------------------------------------------------------
  # Weights follow an asymmetric Dirichlet distribution
  # --------------------------------------------------------------------------
  
  
  # Sample weights from a Dirichlet distribution (with symmetric parameter)
  # alpha <- 1:n / sum(1:n) * n
  #
  alpha <- rep(1, n)
  #
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  # Compute theta_hat for all sampled weights
  theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
  
  # Norms on weights
  w_bar <- alpha / sum(alpha)
  squared_2norm_w_bar <- sum(w_bar^2)
  expected_2norm2 <- mean( apply( samples, 1, function(w) { sum(w^2) } ) )
  expected_3norm3 <- mean( apply( samples, 1, function(w) { sum(w^3) } ) )
  expected_4norm4 <- mean( apply( samples, 1, function(w) { sum(w^4) } ) )
  exact_expected_2norm2 <- sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  exact_expected_3norm3 <- sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) ) )
  exact_expected_4norm4 <- sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 4 ) ) )
  
  expected_2norm2; exact_expected_2norm2 # they should be close
  expected_3norm3; exact_expected_3norm3 # they should be close
  expected_4norm4; exact_expected_4norm4 # they should be close
  
  
  
  # Compute RHS of identity
  # debugonce( M4Exact )
  RHS <- M4Exact( z = z,
                  w_bar = w_bar,
                  exp4norm4 = exact_expected_4norm4,
                  expw3 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) ),
                  expw2 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  RHS  
  sum( (theta_hat - mean(theta_hat))^4 / length(theta_hat) )
  
  
  
    
  
  
  
  
  
  
  
  
  
  
  