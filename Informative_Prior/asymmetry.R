  
  
  ############################################################################
  ### EXACT MOMENTS - INFORMATIVE PRIOR - ASYMMETRY
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     15.05.2022
  # First version:    15.05.2022
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
  
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 18
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1:10)
  df <- 5
  x <- rt( n = n, df = df )
  # z <- x^2
  z <- x
    
  
  
  # --------------------------------------------------------------------------
  # Difference in moments
  # --------------------------------------------------------------------------
  
  
  betaFUN <- function(a, m = 2) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
  n <- 10^2
  set.seed(1:10)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  w_bar <- 1:n / sum(1:n)
  alpha <- w_bar * n  
  
  
  # M1
  m1 <- mean(z)
  m1_w <- sum( w_bar * z )
  
  
  # M2
  m2_b <- M2Exact( z = z, exp2norm2 = 2/n - 1/n^2 )
  m2_bb <- M2Exact( z = z, exp2norm2 = 2/(n + 1) )
  beta_moments <- unlist( lapply( w_bar * n, FUN = betaFUN ) ) 
  m2_bb_w <- M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments) )
  m2_bb
  m2_bb_w
  m2_tilde <- sum( (z - m1)^2 ) / n
  m2_tilde_w <- sum( w_bar * (z - m1_w)^2 )
  m2_tilde_w / m2_tilde; m2_bb_w / m2_bb # same same
  
  
  # M3
  m3_bb <- M3BBFUN( x = z, M2 = 1, scl = FALSE )
  m3_bb_w <- M3BBFUN( x = z, M2 = 1, wghts = w_bar, scl = FALSE )
  m3_bb
  m3_bb_w
  m3_tilde <- sum( (z - m1)^3 ) / n
  m3_tilde_w <- sum( w_bar *  (z - m1_w)^3 )
  m3_tilde
  m3_tilde_w
  m3_tilde_w / m3_tilde; m3_bb_w / m3_bb # same same
  
  
  # M4
  m4_bb <- M4BBFUN( x = z, M2 = 1, scl = FALSE )
  m4_bb_w <- M4BBFUN( x = z, M2 = 1, wghts = w_bar, scl = FALSE )
  m4_bb
  m4_bb_w
  m4_tilde <- sum( (z - m1)^4 ) / n
  m4_tilde_w <- sum( w_bar *  (z - m1_w)^4 )
  m4_tilde
  m4_tilde_w
  (6 * m4_tilde_w + 3 * n * m2_tilde_w^2) / (6 * m4_tilde + 3 * n * m2_tilde^2); m4_bb_w / m4_bb  # same same
  
  
  
  
  # Question: Is E( (w - \bar{w})^m ) the same as for case 0 ?
  # Answer: No.
  
  n_vec <- 3:100
  m_vec <- 2:6
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
  
  for ( m in sq(along = m_vec) ) {
    
    for ( i in seq(along = n_vec) ) {
      
      n <- n_vec[i]
      w_bar <- rep(1/n, n)
      alpha <- w_bar * n
      # sum ( unlist( lapply( alpha, FUN = betaFUN) ) ) - sum(w_bar^m)
      2 / (n + 1) - 1 / n
      # P <- rdirichlet( n = 10^4, alpha = rep(1, n) )
      # mean( apply( P, 1, function(x) { sum( (x - w_bar)^2 ) } ) )
      Alpha <- rdirichlet( n = 10^2, alpha = rep(1, n) ) * n
      
      for ( k in 1:nrow(Alpha) )  {
        
        alpha <- as.numeric(Alpha[k, ] )
        beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
        exact_expected_m_norm_m <- sum(beta_moments)  
        tmp <- exact_expected_m_norm_m - sum((alpha/n)^m)
        tmp / biasCorrection( w_bar = alpha/n, m = m )
        tmp * biasCorrection( w_bar = alpha/n, m = m )
        
      }
    }
  }
  
  
