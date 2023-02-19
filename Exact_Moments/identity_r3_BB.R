  
  
  ############################################################################
  ### IDENTITY FOR r=3, DIRICHLET CASE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.07.2022
  # First version:    20.07.2022
  # --------------------------------------------------------------------------
  
  
  # Description:
  # This script defines an (arbitrarily) real-valued vector z and computes the 
  # RHS and LHS of the identity for the case:
  # - Weights follow a Dirichlet distribution.
  
  
  
  # --------------------------------------------------------------------------
  # Source
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(gtools)
  
  
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Functions
  # --------------------------------------------------------------------------
  
  
  # # --------------------------------------------------------------------------
  # dirichletProductMoments <- function( alpha, type = c("i3", "i2j", "ijk", 
  #                                                      "i4", "i3j", "i2j2", "i2jk", "ijkl") )
  # {
  #   type <- match.arg( type )
  #   ans <- eval( parse( text = paste0("dirichletProductMoments.", type, "(alpha = alpha)" ) ) )
  #   return( ans )
  # }
  # # --------------------------------------------------------------------------
  # dirichletProductMoments.i3 <- function( alpha )
  # {
  #   alpha0 <- sum(alpha)
  #   denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
  #   ans <- (alpha + 2) * (alpha + 1) * alpha / denom
  #   return( ans )
  # }
  # # --------------------------------------------------------------------------
  # dirichletProductMoments.i2j <- function( alpha )
  # {
  #   alpha0 <- sum(alpha)
  #   denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
  #   n <- length(alpha)
  #   ans <- matrix( NA, nrow = length(alpha), ncol = length(alpha) )
  #   for ( i in 1:n ) {
  #     for ( j in 1:n ) {
  #       if( i != j ) {
  #         ans[i, j] <- ( (alpha[i] + 1) * alpha[i] * alpha[j] ) / denom 
  #       }
  #     }
  #   }
  #   return( ans )
  # }
  # # --------------------------------------------------------------------------
  # dirichletProductMoments.ijk <- function( alpha )
  # {
  #   alpha0 <- sum(alpha)
  #   denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
  #   n <- length(alpha)
  #   ans <- array( NA, dim = c(n, n, n) )
  #   for ( i in 1:n ) {
  #     for ( j in 1:n ) {
  #       for ( k in 1:n ) {
  #         if ( i != j && i != k && j != k ) {
  #           ans[i,j,k] <- prod( alpha[c(i, j, k)] ) / denom 
  #         }
  #       }
  #     }
  #   }
  #   return( ans )
  # }
  # dirichletProductMoments( alpha, type = "i3" )
  # dirichletProductMoments( alpha, type = "i2j" )
  # dirichletProductMoments( alpha, type = "ijk" )
  
  
  
  
  # --------------------------------------------------------------------------
  exactEwi2wj <- function( alpha, i, j )
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- ( (alpha[i] + 1) * alpha[i] * alpha[j] ) / denom 
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  exactEwiwjwk <- function(alpha, i, j, k)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- prod( alpha[c(i, j, k)] ) / denom 
    return( ans )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Data and parameters
  # --------------------------------------------------------------------------
  
  # Define real-valued vector z
  n <- 4
  set.seed(1111)
  z <- rt(n = n, df = 3)^2
  
  
  
  # Number of simulations (for the Dirichlet sampling case)
  n_sim <- 10^6
  
  
 
  

  
  # --------------------------------------------------------------------------
  # Weights follow a flat Dirichlet distribution
  # Conclusion: LHS = RHS
  # --------------------------------------------------------------------------
  
  # Sample weights from a Dirichlet distribution
  lambda <- 1
  alpha <- rep(1, n) * lambda 
  w_bar <- alpha / sum(alpha)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar
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
  sum( (theta_hat - mean(theta_hat))^3 ) / length(theta_hat) # close enough
  
  
  
  
  
  
  
  
  # Compute LHS using E(w_i^3), E(w_i^2 w_j) and E(w_i w_j, w_k)
  alpha0 <- sum(alpha)
  
  wi3_bar <- unlist( lapply( alpha, FUN = betaMoments, a0 = sum(alpha), m = 3 ) )
  ( (alpha[1] + 2) * (alpha[1] + 1) * alpha[1] ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  wi3_bar # same same
  
  w1w2 <- mean( apply( samples, 1, function(w) { prod(w[1:2] ) } ) )
  ( alpha[1] * alpha[2] ) / (alpha0 * (alpha0 + 1) )
  w1w2 # same same
  
  w1sqw2 <- mean( apply( samples, 1, function(w) { w[1]^2 * w[2] } ) )
  ( (alpha[1] + 1) * alpha[1] * alpha[2] ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  w1sqw2 # same same
  
  
 
  
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
    ans <- apply( A, 2, function(a) { ( (a[1] + 1) * a[1] * a[2] ) / denom } )
    return( ans )
  }
  
  exactEwi2wj( alpha = alpha )
  wi2wj_bar
  
  wi2wj_mat <- matrix(0, n, n)
  wi2wj_mat[upper.tri(wi2wj_mat)] <- wi2wj_bar
  wi2wj_mat <- t(wi2wj_mat)
  wi2wj_mat[upper.tri(wi2wj_mat)] <- wi2wj_bar
  wi2wj_mat <- t(wi2wj_mat)
  wi2wj_mat
  
  
  
  
  
  
  
  
  # Compute LHS
  wi2wj_mat <- matrix(0, n, n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        wi2wj_mat[i, j] <- exactEwi2wj( alpha = alpha, i = i, j = j )
      }
    }
  }
  wiwjwk_array <- array(0, dim = c(n, n, n))
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        wiwjwk_array[i, j, k] <- exactEwiwjwk( alpha = alpha, i = i, j = j, k = k )
      }
    }
  }
  alpha0 <- sum(alpha)
  wi3_vec <- ( (alpha + 2) * (alpha + 1) * alpha ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  tmp1 <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        tmp1 <- tmp1 + wi2wj_mat[i, j] * b[i]^2  * b[j] * 3
      }
    }
  }
  tmp2 <- 0
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        tmp2 <- tmp2 + wiwjwk_array[i, j, k] * prod(b[c(i, j, k)]) * 6
      }
    }
  }
  LHS <- sum(wi3_vec * b^3) + tmp1 + tmp2
  LHS
  RHS # same same
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Weights follow a symmetric Dirichlet distribution, \alpha_0 <> n
  # Conclusion: LHS = RHS
  # --------------------------------------------------------------------------
  
  # Sample weights from a Dirichlet distribution (with symmetric parameter)
  lambda <- 0.1
  alpha <- rep(1, n) * lambda
  w_bar <- alpha / sum(alpha)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  cbind( w_bar, apply( samples, 2, mean) )
  
  
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
  sum( (theta_hat - mean(theta_hat))^3 ) / length(theta_hat)
  
  
  
  
  
  
  
  alpha0 <- sum(alpha)
  wi3_vec <- ( (alpha + 2) * (alpha + 1) * alpha ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  wi2wj_mat <- matrix(0, n, n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        wi2wj_mat[i, j] <- exactEwi2wj( alpha = alpha, i = i, j = j )
      }
    }
  }
  wiwjwk_array <- array(0, dim = c(n, n, n))
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        wiwjwk_array[i, j, k] <- exactEwiwjwk( alpha = alpha, i = i, j = j, k = k )
      }
    }
  }
  tmp1 <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        tmp1 <- tmp1 + wi2wj_mat[i, j] * b[i]^2  * b[j] * 3
      }
    }
  }
  tmp2 <- 0
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        tmp2 <- tmp2 + wiwjwk_array[i, j, k] * prod(b[c(i, j, k)]) * 6
      }
    }
  }
  LHS <- sum(wi3_vec * b^3) + tmp1 + tmp2
  LHS
  RHS # same same
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Weights follow an asymmetric Dirichlet distribution, alpha sums to n
  # Conclusion: LHS = RHS
  # --------------------------------------------------------------------------
  
  
  # Sample weights from a Dirichlet distribution (with asymmetric parameter)
  alpha <- 1:n / sum(1:n) * n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  w_bar <- alpha / sum(alpha)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar  # important step!
  cbind( w_bar, apply( samples, 2, mean) )
  
  # Compute theta_hat for all sampled weights
  theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
  
  # Norms on weights
  expected_2norm2 <- mean( apply( samples, 1, function(w) { sum(w^2) } ) )
  expected_3norm3 <- mean( apply( samples, 1, function(w) { sum(w^3) } ) )
  exact_expected_2norm2 <-  sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  exact_expected_3norm3 <- sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) ) )

  expected_2norm2; exact_expected_2norm2  # should be similar
  expected_3norm3; exact_expected_3norm3  # should be similar
 
  
  
  # Compute RHS of identity
  RHS <- M3Exact( z = z,
                  w_bar = w_bar,
                  exp3norm3 = exact_expected_3norm3,
                  # exp2norm2 = exact_expected_2norm2,
                  expw2 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  RHS 
  sum( (theta_hat - mean(theta_hat))^3 ) / length(theta_hat) # similar
  M3BBFUN( x = z, M2 = Inf, wghts = w_bar, scl = FALSE ) # same
  2 / ((alpha0 + 2) * (alpha0 + 1)) * sum( alpha / alpha0 * b^3 ) # same
  
  
  
  
  # Compute LHS
  alpha0 <- sum(alpha)
  wi3_vec <- ( (alpha + 2) * (alpha + 1) * alpha ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  wi2wj_mat <- matrix(0, n, n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        wi2wj_mat[i, j] <- exactEwi2wj( alpha = alpha, i = i, j = j )
      }
    }
  }
  wiwjwk_array <- array(0, dim = c(n, n, n))
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        wiwjwk_array[i, j, k] <- exactEwiwjwk( alpha = alpha, i = i, j = j, k = k )
      }
    }
  }
  tmp1 <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        tmp1 <- tmp1 + wi2wj_mat[i, j] * b[i]^2  * b[j] * 3
      }
    }
  }
  tmp2 <- 0
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        tmp2 <- tmp2 + wiwjwk_array[i, j, k] * prod(b[c(i, j, k)]) * 6
      }
    }
  }
  LHS <- sum(wi3_vec * b^3) + tmp1 + tmp2
  LHS
  RHS # same same
  
  
  
    
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Weights follow an asymmetric Dirichlet distribution, \alpha_0 <> n
  # Conclusion: LHS = RHS
  # --------------------------------------------------------------------------
  
  
  # Sample weights from a Dirichlet distribution (with asymmetric parameter)
  alpha <- 1:n / sum(1:n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  w_bar <- alpha / sum(alpha)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar  # important step!
  cbind( w_bar, apply( samples, 2, mean) )
  
  # Compute theta_hat for all sampled weights
  theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
  
  # Norms on weights
  expected_2norm2 <- mean( apply( samples, 1, function(w) { sum(w^2) } ) )
  expected_3norm3 <- mean( apply( samples, 1, function(w) { sum(w^3) } ) )
  exact_expected_2norm2 <-  sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  exact_expected_3norm3 <- sum( unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) ) )
  
  expected_2norm2; exact_expected_2norm2  # should be similar
  expected_3norm3; exact_expected_3norm3  # should be similar
  
  
  
  # Compute RHS of identity
  RHS <- M3Exact( z = z,
                  w_bar = w_bar,
                  exp3norm3 = exact_expected_3norm3,
                  # exp2norm2 = exact_expected_2norm2,
                  expw2 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  RHS 
  sum( (theta_hat - mean(theta_hat))^3 ) / length(theta_hat) # similar
  M3BBFUN( x = z, M2 = Inf, wghts = w_bar, scl = FALSE ) # same
  2 / ((alpha0 + 2) * (alpha0 + 1)) * sum( alpha / alpha0 * b^3 ) # same
  
  
  
  
  # Compute LHS
  alpha0 <- sum(alpha)
  wi3_vec <- ( (alpha + 2) * (alpha + 1) * alpha ) / ( (alpha0 + 2) * (alpha0 + 1) * alpha0 )
  wi2wj_mat <- matrix(0, n, n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        wi2wj_mat[i, j] <- exactEwi2wj( alpha = alpha, i = i, j = j )
      }
    }
  }
  wiwjwk_array <- array(0, dim = c(n, n, n))
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        wiwjwk_array[i, j, k] <- exactEwiwjwk( alpha = alpha, i = i, j = j, k = k )
      }
    }
  }
  tmp1 <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        tmp1 <- tmp1 + wi2wj_mat[i, j] * b[i]^2  * b[j] * 3
      }
    }
  }
  tmp2 <- 0
  for ( i in 1:(n-2) ) {
    for ( j in (i+1):(n-1) ) {
      for ( k in (j+1):n ) {
        tmp2 <- tmp2 + wiwjwk_array[i, j, k] * prod(b[c(i, j, k)]) * 6
      }
    }
  }
  LHS <- sum(wi3_vec * b^3) + tmp1 + tmp2
  LHS
  RHS # same same
  
  
  
  
  
  
  