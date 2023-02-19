  
  
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
  # RHS and LHS of the identity for the case where Weights follow a Dirichlet distribution. 
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Source
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(gtools)
  require(Umoments)
  
  
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
  exactEwi4 <- function(alpha, i, j)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 3) * (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- ( (alpha[i] + 3) * (alpha[i] + 2) * (alpha[i] + 1) * alpha[i] ) / denom 
    return( ans )
  }
  # --------------------------------------------------------------------------
  exactEwi3wj <- function(alpha, i, j)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 3) * (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- ( (alpha[i] + 2) * (alpha[i] + 1) * alpha[i] * alpha[j] ) / denom 
    return( ans )
  }
  # --------------------------------------------------------------------------
  exactEwi2wj2 <- function(alpha, i, j)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 3) * (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- ( (alpha[i] + 1) * alpha[i] * (alpha[j] + 1) * alpha[j] ) / denom 
    return( ans )
  }
  # --------------------------------------------------------------------------
  exactEwi2wjwk <- function(alpha, i, j, k)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 3) * (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- ( (alpha[i] + 1) * prod(alpha[c(i,j,k)]) ) / denom 
    return( ans )
  }
  # --------------------------------------------------------------------------
  exactEwiwjwkwl <- function(alpha, i, j, k, l)
  {
    alpha0 <- sum(alpha)
    denom <- (alpha0 + 3) * (alpha0 + 2) * (alpha0 + 1) * alpha0
    ans <- prod(alpha[c(i,j,k,l)]) / denom 
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  # Data and parameters
  # --------------------------------------------------------------------------
  
  # Define real-valued vector z
  n <- 10
  set.seed(1111)
  z <- rt(n = n, df = 3)^2
  
  # Number of simulations (for the Dirichlet sampling case)
  n_sim <- 10^6
  
  
 
  

  
  # --------------------------------------------------------------------------
  # Weights follow a flat Dirichlet distribution
  # --------------------------------------------------------------------------
  
  # Sample weights from a Dirichlet distribution
  alpha <- rep(1, n) 
  w_bar <- alpha / sum(alpha)
  z_bar <- sum( w_bar * z )
  b <- z - z_bar
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  # Compute theta_hat for all sampled weights
  theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
  
  # Norms on weights
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
  RHS <- M4Exact( z = z,
                  w_bar = w_bar,
                  exp4norm4 = exact_expected_4norm4,
                  expw3 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) ),
                  expw2 = unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) ) )
  RHS
  sum( (theta_hat - mean(theta_hat))^4 ) / length(theta_hat)
  
  
  
  
  scl4_w <- biasCorrection( N = n, m = 4 )
  scl4 <- biasCorrection( N = n, m = 4, unweighted = TRUE )
  scl4_w[1:2] / scl4_w[3]
  scl4[1:2] / scl4[3]

  
  
  
  
  
  # Compute LHS using E(w_i^3), E(w_i^3 w_j), E(w_i^2 w_j w_k) and E(w_i w_j w_k w_l)
 
  wi4_vec <- unlist( lapply( 1:n, FUN = exactEwi4, alpha = alpha ) )
  wi3wj_mat <- matrix(0, n, n)
  wi3wj_mat_est <- wi3wj_mat
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        wi3wj_mat[i, j] <- exactEwi3wj( alpha = alpha, i = i, j = j )
        wi3wj_mat_est[i, j] <- mean( apply( samples, 1, function(x) { x[i]^3 * x[j] } ) )
      }
    }
  }
  wi3wj_mat
  wi3wj_mat_est
  
  wi2wj2_mat <- matrix(0, n, n)
  wi2wj2_mat_est <- wi2wj2_mat
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      wi2wj2_mat[i, j] <- exactEwi2wj2( alpha = alpha, i = i, j = j )
      wi2wj2_mat_est[i, j] <- mean( apply( samples, 1, function(x) { x[i]^2 * x[j]^2 } ) )
    }
  }
  wi2wj2_mat
  wi2wj2_mat_est
  
  wi2wjwk_array <- array(0, dim = c(n, n, n))
  wi2wjwk_array_est <- wi2wjwk_array
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      for ( k in 1:n ) {
        if ( i != j && j != k && i != k ) {
          wi2wjwk_array[i, j, k] <- exactEwi2wjwk( alpha = alpha, i = i, j = j, k = k )
          wi2wjwk_array_est[i, j, k] <- mean( apply( samples, 1, function(x) { x[i]^2 * x[j] * x[k] } ) )
        }
      }
    }
  }
  wi2wjwk_array
  wi2wjwk_array_est
  
  wiwjwkwl_array <- array(0, dim = c(n, n, n, n))
  wiwjwkwl_array_est <- wiwjwkwl_array
  for ( i in 1:(n-3) ) {
    for ( j in (i+1):(n-2) ) {
      for ( k in (j+1):(n-1) ) {
        for ( l in (k+1):n ) {
          wiwjwkwl_array[i,j,k,l] <- exactEwiwjwkwl( alpha = alpha, i = i, j = j, k = k, l = l )
          wiwjwkwl_array_est[i,j,k,l] <- mean( apply( samples, 1, function(x) { prod(x) } ) )
        }
      }
    }
  }
  
  tmp1 <- 0
  tmp1_est <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i != j ) {
        tmp1 <- tmp1 + wi3wj_mat[i, j] * b[i]^3 * b[j] * 4
        tmp1_est <- tmp1_est + wi3wj_mat_est[i, j] * b[i]^3 * b[j] * 4
      }
    }
  }
  tmp2 <- 0
  tmp2_est <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp2 <- tmp2 + wi2wj2_mat[i, j] * b[i]^2 * b[j]^2 * 6
      tmp2_est <- tmp2_est + wi2wj2_mat_est[i, j] * b[i]^2 * b[j]^2 * 6
    }
  }
  tmp3 <- 0
  tmp3_est <- 0
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      for ( k in 1:n ) {
        if ( i != j && j != k && i != k ) {
          tmp3 <- tmp3 + wi2wjwk_array[i, j, k] * b[i]^2 * b[j] * b[k] * 6
          tmp3_est <- tmp3_est + wi2wjwk_array_est[i, j, k] * b[i]^2 * b[j] * b[k] * 6
        }
      }
    }
  }
  tmp4 <- 0
  tmp4_est <- 0
  for ( i in 1:(n-3) ) {
    for ( j in (i+1):(n-2) ) {
      for ( k in (j+1):(n-1) ) {
        for ( l in (k+1):n ) {
          tmp4 <- tmp4 + wiwjwkwl_array[i,j,k,l] * prod(b[c(i,j,k,l)]) * 24
          tmp4_est <- tmp4_est + wiwjwkwl_array_est[i,j,k,l] * prod(b[c(i,j,k,l)]) * 24
        }
      }
    }
  }
  LHS <- sum(wi4_vec * b^4) + tmp1 + tmp2 + tmp3 + tmp4
  LHS
  sum(wi4_vec * b^4) + tmp1_est + tmp2_est + tmp3_est + tmp4_est # similar enough - not so for n = 10
  M4BBFUN( x = z, M2 = Inf, scl = FALSE ) # same same
  sum( (theta_hat - mean(theta_hat))^4 ) / length(theta_hat) # similar enough
  RHS 
  
 
  
  
  
  
  
  alpha <- rep(1, n) 
  w_bar <- alpha / sum(alpha)
  expw2 <- unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 2 ) )
  expw3 <- unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 3 ) )
  expw4 <- unlist( lapply( alpha, betaMoments, a0 = sum(alpha), m = 4 ) )
  exact_expected_4norm4 <- sum( expw4 )
  scl4 <- biasCorrection( w_bar = w_bar, m = 4 )
  scl4_alt <- biasCorrection( N = length(z), m = 4, unweighted = TRUE )
  D <- ( exact_expected_4norm4 - 3 * sum(w_bar^4) - 4 * sum(w_bar * expw3) + 6 * sum(w_bar^2 * expw2) )
  D2 <- sum(expw2) - sum(w_bar^2)
  D
  mean( apply( samples, 1, function(x) { sum((x - w_bar)^4) } ) )
  
  N <- 10^2
  rhs <- rep(NA, N)
  # rhs_est <- rhs
  lhs <-  rep(NA, N)
  m4 <- rep(NA, N)
  m4_alt <- m4
  m4_uM <- m4
  m4_ml_vec <- m4
  m2_ml_vec <- m4
  a_vec <- m4
  b_vec <- m4
  seed_vec <- unlist( lapply( 1:N, FUN = function(i) { runif(1, 1, 10^7) } ) )
  
  for ( j in seq(along = lhs) ) {
    
    set.seed( seed_vec[j] )
    # z <- rt(n = n, df = 4)
    # z <- runif(n = n, 1, 10)
    z <- rnorm(n = n)
    
    # samples <- rdirichlet( n = n_sim, alpha = alpha )
    # theta_hat <- apply( samples, 1, function(w) { sum(w * z) } )
    # rhs_est[j] <- sum( (theta_hat - mean(theta_hat))^4 ) / length(theta_hat)
    lhs[j] <-  M4BBFUN( x = z, M2 = Inf, scl = FALSE )
    rhs[j] <- M4Exact( z = z,
                       w_bar = w_bar,
                       exp4norm4 = exact_expected_4norm4,
                       expw3 = expw3,
                       expw2 = expw2 )
    z_bar <- sum( w_bar * z )
    m2_ml <- sum( w_bar * (z - z_bar)^2 )
    m4_ml <- sum( w_bar * (z - z_bar)^4 )
    # m4[j] <- (m4_ml * scl4["scl4a"] - m2_ml^2 * scl4["scl4b"]) / scl4["scl4denom"]
    m4[j] <- (m4_ml * scl4["scl4a"] + (m2_ml)^2 * scl4["scl4b"]) / n^3  #(scl4["scl4denom"] * n)
    m4_alt[j] <- (m4_ml * scl4_alt["scl4a"] + m2_ml^2 * scl4_alt["scl4b"]) / scl4_alt["scl4denom"]
    a_vec[j] <- m4_ml * scl4_alt["scl4a"] / scl4_alt["scl4denom"]
    b_vec[j] <- m2_ml^2 * scl4_alt["scl4b"] / scl4_alt["scl4denom"]
    m4_uM[j] <- uM4( m2 = m2_ml, m4 = m4_ml, n = n )
    m4_ml_vec[j] <- m4_ml
    m2_ml_vec[j] <- m2_ml
    
  }
  
  # plot( x = rhs_est, y = rhs, ylim = range(c(rhs, rhs_est)), xlim = range(c(rhs, rhs_est)) )
  # points( x = rhs_est, y = rhs_est, pch = 19, col = 2 )
  # cbind( rhs_est, rhs )
  
  
  plot( m4_uM )
  which( m4_uM < 0 )
  
  
  
  
  
  plot( x = lhs, y = rhs, ylim = range(c(lhs, rhs)), xlim = range(c(lhs, rhs)) )
  points( x = lhs, y = lhs, pch = 19, col = 2 )
  points( x = lhs, y = (m4 * D), pch = 19, col = 3 )
  points( x = lhs, y = m4_alt * D, pch = 19, col = 4 )
  points( x = lhs, y = (m4_uM * D), pch = 19, col = 5 )
  points( x = lhs, y = (m4_ml_vec * D), pch = 19, col = 6 )
  points( x = lhs, y = (6*m4_ml_vec + 3*n*m2_ml_vec^2)/((n+3)*(n+2)*(n+1)), pch = 19, col = 7 )
  points( x = lhs, y = a_vec * D + b_vec * D2^2, pch = 19, col = 8 )
  
  
  
  fac <- (max(lhs) - (min(lhs) - min(m4*D))) / max(m4*D)
  points( x = lhs, y = m4*D*fac, col = 8 )
  
  
  tmp <- a_vec * D - b_vec * (-1) * D2^2
  fac2 <- (max(lhs) - (min(lhs) - min(tmp))) / max(tmp)
  points( x = lhs, y = tmp * fac2, col = 9, pch = 19 )
  
  
  
  
  
  
  
  cbind( lhs, m4 * D, m4 * D * fac, m4_uM * D, m4_alt * D, m4_ml * D )
  
  min(lhs); min(m4*D*fac)
  max(lhs); max(m4*D*fac)
  
  
  
  m4_should <- lhs / D
  
  fac <- lhs / (m4 * D)
  range(fac)
  
  
  
  
  m4
  m4_alt
  plot(m4_alt)
  
  
  plot( x = m4, y = m4_alt )
  
  
  
  lhs / m4
  range( lhs / m4 )
  D
  
  
   
  

  
  
  #####################
  
  # install.packages("Umoments")
  require(Umoments)  
  
  # z <- runif(n = n) * 10
  
  uM(z, 4)
  
  b <- z - mean(z)
  mean(b^2) * (n / (n-1))  
  mean(b^3) * n^2 / ((n-1) * (n-2))
    
  
  scl4 <- biasCorrection( N = length(z), m = 4, unweighted = TRUE )
  m2_ml <- mean(b^2)
  m3_ml <- mean(b^3)
  m4_ml <- mean(b^4)
  m4 <- (m4_ml * scl4["scl4a"] + m2_ml^2 * scl4["scl4b"]) / scl4["scl4denom"]
  m4
  
  
  # uM3 <- function(m3, n) 
  # {
  #   m3 * n^2 / ((n - 1) * (n - 2))
  # }
  # 
  # uM4 <- function(m2, m4, n) 
  # {
  #   -3 * m2^2 * (2 * n - 3) * n/((n - 1) * (n - 2) * (n - 3)) + 
  #     (n^2 - 2 * n + 3) * m4 * n/((n - 1) * (n - 2) * (n - 3))
  # }
  
  
  Umoments::uM3( m3 = m3_ml, n = n )
  Umoments::uM4( m2 = m2_ml, m4 = m4_ml, n = n )
  
  
  
  
 
  
  ##########################################################
  
  N <- n
  denom <- (N - 1) * (N - 2) * (N - 3)
  a1 <- N * (N^2 - 2*N + 3) / denom
  a2 <- 3*N * (3 - 2*N) / denom
  # a2 <- -3*N * (2*N - 3) / denom
  # a2 <- -3*N * (3 - 2*N) / denom
  b1 <- 6 / ((N+1)*(N+2)*(N+3))
  b2 <- 3*N / ((N+1)*(N+2)*(N+3))

  b1 / a1  
  b2 / a2

  a1 * (b1/a1); b1
  a2 * (b2/a2); b2
  
   
  
  
  # set.seed(1111)
  # z <- runif(n = n, 1, 10)
  w_bar <- rep(1/n, n)
  m2_ml <- mean( (z - sum(w_bar * z))^2 )
  m4_ml <- mean( (z - sum(w_bar * z))^4 )
  (a1 * m4_ml + a2 * m2_ml^2) * D + (a1 * m4_ml + a2 * m2_ml^2) * (sum(expw2) - sum(w_bar^2))

  LHS
  
  
    
  
  
  
  
  # --------------------------------------------------------------------------
  # CHECK CONSISTENCY AND BIAS OF UNBIASED MOMENT ESTIMATORS
  # --------------------------------------------------------------------------
  
  v <- 5
  um4 <- numeric(10^2)
  for ( j in seq(along = um4 ) ) {
    r <- rnorm( n = 4 )
    # r <- rt( n = 4, df = v )
    m2_ml <- mean( (r - mean(r))^2 )
    m4_ml <- mean( (r - mean(r))^4 )
    um4[j] <- Umoments::uM4( m2 = m2_ml, m4 = m4_ml, n = length(r) )
  }
  
  range(um4); mean(um4)
  
  
  v / (v - 2); m2_ml
  
  range(um4)
  mean(um4)
  (6 / (v - 4) + 3) * (v / (v - 2))^2
  
  
  
  
  
  
  
  
  