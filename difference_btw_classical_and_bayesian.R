  
  
  ############################################################################
  ### EXACT MOMENTS - DIFFERENCE BETWEEN CLASSICAL AND BAYESIAN BOOTSTRAP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.05.2022
  # First version:    13.05.2022
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
  
  
  
  # Sections:
  # Difference in moments, loop over increasing n
  # Mapping between m out of n bootstrap and lambda parameter in symmetric BB 
  
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  # n <- 7
  # n_boot <- 50
  # B <- 10^4 + 1
  # set.seed(1111)
  # df <- 5
  # x <- rt( n = n, df = df )
  # z <- x^2

  
  
  # Bertsimas Sturt 2020 example
  
  n <- 81
  B <- 10^5
  n_boot <- 10^1
  p_th <- 0.025
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, n - length(z)))
  z_demeaned <- scale(z, TRUE, FALSE)
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
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
  
  
  
  
  # --------------------------------------------------------------------------
  # Difference in moments, loop over increasing n
  # --------------------------------------------------------------------------
  
  # M2
  #// Me: difference converges if z has finite variance
 
  n_vec <- 3:10^3
  set.seed(1111)
  z_vec <-  rt( n = max(n_vec), df = 4 )^2
  # z_vec <- n_vec^3
  M2 <- matrix( NA, nrow = length(n_vec), ncol = 6, 
                dimnames = list(NULL, c("B", "BB", "m2", "m2_tilde", "EB", "EBB")) )
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    z <- z_vec[1:n]
    m2 <- var(z)
    m2_tilde <- sum( (z - mean(z))^2 ) / n
    m2_b <- M2Exact( z = z, exp2norm2 = 2 / n - 1 / n^2 )
    m2_bb <- M2Exact( z = z, exp2norm2 = 2 / (n + 1) )
    M2[i, ] <- c(m2_b, m2_bb, m2, m2_tilde, 2/n - 1/n^2, 2 / (n + 1))
  }
  
  a <- M2[ ,"m2"] * (M2[ ,"EB"] - M2[ ,"EBB"])
  b <- M2[ ,"m2_tilde"] * (1/n_vec - 1/(n_vec+1))
  plot( a - b )                                         # same same
  plot( b - M2[ ,"m2_tilde"] * 1/(n_vec^2 + n_vec) )    # same same
  
  plot( a )
  
  
  plot( M2[ ,"B"] - M2[ ,"BB"] )
  abline( h = 0 )
  
  
  plot( M2[ ,"EBB"] - M2[ ,"EB"] )
  
  plot( 1 / (n_vec + 1) - 1 / n_vec )
  plot( (1 / (n_vec + 1) - 1 / n_vec) * M2[ ,"m2_tilde"] )
  
  
 

    
  
  
  # M3
  
  M3 <- matrix( NA, nrow = length(n_vec), ncol = 5, 
                dimnames = list(NULL, c("B", "BB", "m3_tilde", "EB", "EBB")) )
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    z <- z_vec[1:n]
    m3_tilde <- sum( (z - mean(z))^3 ) / n
    m3_b <- M3BFUN( x = z, M2 = NULL, scl = FALSE )
    m3_bb <- M3BBFUN( x = z, M2 = 1, scl = FALSE )
    M3[i, ] <- c(m3_b, m3_bb, m3_tilde, 2/n - 1/n^2, 2 / (n + 1))
  }
  
  head(M3)
  
  plot( M3[ ,"BB"] - M3[ ,"B"] )
  abline( h = 0 )
  
  plot( (M3[ ,2] - M3[ ,1]) / M3[ ,3] )
  abline( h = 0 )
  
  plot( M2[ ,"EBB"] - M2[ ,"EB"] )
  
  plot( 1 / (n_vec + 1) - 1 / n_vec )
  plot( (1 / (n_vec + 1) - 1 / n_vec) * M2[ ,"m2_tilde"] )
  
  
  
  # M4
  
  M4 <- matrix( NA, nrow = length(n_vec), ncol = 5, 
                dimnames = list(NULL, c("B", "BB", "m3_tilde", "EB", "EBB")) )
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    z <- z_vec[1:n]
    m4_tilde <- sum( (z - mean(z))^4 ) / n
    m4_b <- M4BFUN( x = z, M2 = NULL, scl = FALSE )
    # debugonce( M4BBFUN )
    m4_bb <- M4BBFUN( x = z, M2 = 1, scl = FALSE )
    M4[i, ] <- c(m4_b, m4_bb, m4_tilde, 2/n - 1/n^2, 2 / (n + 1))
  }
  
  plot( M4[ ,"BB"] - M4[ ,"B"] )
  abline( h = 0 )  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Mapping between m out of n bootstrap and lambda parameter in symmetric BB 
  # ( Dir(alpha * lambda) )
  # --------------------------------------------------------------------------
  
  
  mapFUN <- function(n, m) { (1-m)/(1-n) + (m-1)/(n*(1-n)) }
  n_vec <- 1:100
  m_vec <- 1:100
  lambda_mat <- matrix( NA, nrow = length(n_vec), ncol = length(m_vec) )
  e2n2_B <- lambda_mat
  e2n2_BB <- lambda_mat
  for ( i in seq(along = n_vec) ) {
    for ( j in seq(along = m_vec) ) {
      n <- n_vec[i]
      m <- m_vec[j]
      lambda_mat[i, j] <- mapFUN( n = n, m = m )
      e2n2_B[i, j] <- (n + m - 1) / (n * m)
      e2n2_BB[i, j] <- (1 + lambda_mat[i, j]) / (1 + lambda_mat[i, j] * n)
    }
  }
  e2n2_B - e2n2_BB
  
  
  
  # m = n:
  
  mapFUN <- function(n) { (n^2 - 2*n + 1) / (n^2 - n) } 
  n_vec <- 1:100
  lambda_vec <- rep(NA, length(n_vec))
  e2n2_B <- lambda_vec
  e2n2_BB <- lambda_vec
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    lambda_vec[i] <- mapFUN( n = n )
    e2n2_B[i] <- (2*n - 1) / n^2
    e2n2_BB[i] <- (1 + lambda_vec[i]) / (1 + lambda_vec[i] * n)
  }
  e2n2_B - e2n2_BB
  
  a <- (2*n_vec - 1) / n_vec^2
  b <- (1 + lambda_vec) / (1 + lambda_vec * n_vec)
  cbind( a, b )
  plot( a, b )
  
  
  
  
  # --------------------------------------------------------------------------
  # Mapping between BB and m out of n bootstrap to get the same std. dev.
  # --------------------------------------------------------------------------
  
  n_vec <- 3:100
  m_vec <- 3:100
  
  lambda_mat <- matrix( NA, nrow = length(n_vec), ncol = length(m_vec) )
  lambda_mat_approx <- lambda_mat
  for ( i in seq(along = n_vec) ) {
    for ( j in seq(along = m_vec) ) {
      lambda_mat[i, j] <- mn2Lambda( n = n_vec[i], m = m_vec[j] )
      lambda_mat_approx[i, j] <- (m_vec[j] - 1) / n_vec[i]
    }
  }
  headleft( lambda_mat )
  headleft( lambda_mat_approx )
  
  
  range( lambda_mat - lambda_mat_approx )
  
  
  
  
  