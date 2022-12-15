  

  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
    
  n <- 2
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1:10)
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
  
  
  ewi <- 2/n^2 - 1/n^3
  b <- z - mean(z)
  ewi * b[1]^2 + ewi * b[2]^2 + 2 * w_bar[1] * w_bar[2] * b[1] * b[2]
  (ewi * b[1]^2 + ewi * b[2]^2 + 2 * w_bar[1] * w_bar[2] * b[1] * b[2]) * 2
  

  ( sum( ewi * b^2 ) + 
      2 * prod( b[c(1, 2)] ) * prod( w_bar[c(1, 2)] ) + 
      2 * prod( b[c(1, 3)] ) * prod( w_bar[c(1, 3)] ) +
      2 * prod( b[c(3, 2)] ) * prod( w_bar[c(3, 2)] ) )
      
  m2 * 1/n  
  
  
  
  
  
  
  ####################################################
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  
  tic <- Sys.time()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  mean( apply( boot_mat, 2, var ) )

  
  
  
  M2BFUN <- function( z, wghts = NULL )
  {
    N <- length(z)
    if ( is.null(wghts) ) {
      wghts <- rep(1/N, N)
    }
    z_bar <- sum( wghts * z )
    
    S <- cov.wt( x = matrix(z, ncol = 1), 
                 wt = wghts, 
                 cor = FALSE, 
                 center = z_bar,
                 method = "ML" )
    
    ans <- as.numeric(S$cov) * 1 / N
    
    return( ans )
  }
  
  M2BBFUN <- function( z, wghts = NULL )
  {
    N <- length(z)
    if ( is.null(wghts) ) {
      wghts <- rep(1/N, N)
    }
    z_bar <- sum( wghts * z )
    
    S <- cov.wt( x = matrix(z, ncol = 1), 
                 wt = wghts, 
                 cor = FALSE, 
                 center = z_bar,
                 method = "ML" )
    
    ans <- as.numeric(S$cov) * 1 / (N + 1)
    
    return( ans )
  }
  
  M2BFUN( z = z )
  M2BBFUN( z = z )
  
  
  z_bar <- mean(z)
  S <- var(z)
  S * ( 2/n - 1/n^2 - 1/n ) 
  
  
  
  
  # -------------------------
  # BB
  
  # wghts <- 1:n / sum(1:n)
  wghts <- rep(1/n, n)
  alpha <- wghts * n
  m <- 2
  beta_mom <- unlist( lapply( 1:n, FUN = function(i) { betaMoments( a = alpha[i], b = sum(alpha) - alpha[i], m = m ) } ) )
  sum(beta_mom)    
  2 / (n + 1)
  
  z_bar <- sum(wghts * z)
  S <- cov.wt( x = matrix(z, ncol = 1), 
               wt = wghts, 
               cor = FALSE, 
               center = z_bar,
               method = "unbiased" )
  S$cov * ( sum(beta_mom) - sum(wghts^2) )
  1 / (1 - sum(wghts^2)) * sum( wghts * (z - z_bar)^2 ) * ( sum(beta_mom) - sum(wghts^2) )
  sum( wghts * (z - z_bar)^2 ) / (n+1)  

  
  beta_mom * (z - z_bar)^2
  sum( beta_mom * (z - z_bar)^2 )
  
  b <- z - z_bar
  sum( beta_mom * b^2 ) + 
    2 * prod( b[c(1, 2)] ) * prod( w_bar[c(1, 2)] ) + 
    2 * prod( b[c(1, 3)] ) * prod( w_bar[c(1, 3)] ) +
    2 * prod( b[c(3, 2)] ) * prod( w_bar[c(3, 2)] )
  
  
  
  sum( beta_mom * z^2 ) - sum( wghts^2 * z^2 )
  
  
  
  
  
  
  
  
  
  
  
    