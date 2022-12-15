  
  
  ############################################################################
  ### DISTRIBUTION OF LARGEST EIGENVALUE - SYNTHETIC DATA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     08.08.2022
  # First version:    08.08.2022
  # --------------------------------------------------------------------------
  
  
  # Replica of El Karoui and Purdom 2016:
  # The bootstrap, covariance matricces and PCA in moderate and high-dimensions
  
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(rgl)
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
  
  N <- 10^5
  n <- 10^3
  p <- 10^3 / 2
  
  X <- rmvnorm( nObs = N, nAssets = p, asTS = FALSE )
  covmat <- cov(X)
  E <- eigen(covmat)
  E$values[1:10]
  
  
  lX <- list()
  for ( i in 1:(N/n) ) {
    idx <- i*n
    lX[[i]] <- X[(idx-n+1):idx, ]
  }
  lCov <- lapply( lX, cov )
  lE <- lapply( lCov, FUN = eigen )
  lambda_1 <- unlist( lapply( lE, FUN = function(x) { x$values[1] } ) )
  
  plot( density( lambda_1 ), xlim = c(0, max(lambda_1)) ) 
  abline( v = 1 )
    
  
  
  
  i <- 10
  E <- eigen(lCov[[i]])
  evec <- E$vectors
  
  S <- scale(lX[[i]], TRUE, FALSE) %*% evec
  headleft(cov(S))
  
  z <- S[ ,1]^2
  M2Exact( z = z, exp2norm2 = 2 / (n + 1) )
  var(z) / n
  
  samples <- rdirichlet( n = 10^4, alpha = rep(1, n) )  
  theta <- apply( samples, 1, function(w) { sum(w * z) } )
  var(theta)
  mean(theta)
  
  
  
  plot( density( lambda_1 ), xlim = c(0, max(c(lambda_1, theta))) ) 
  lines( density(theta), col = 2 )
  abline( v = 1 )
  abline( v = mean(lambda_1) )
  abline( v = mean(theta), col = 2 )
  
    
  
  
  
  FUN <- function(i)
  {
    S <- scale(lX[[i]], TRUE, FALSE) %*% lE[[i]]$vectors
    z <- S[ ,1]^2
    m1 <- mean(z)
    m2 <- M2Exact( z = z, exp2norm2 = 2 / (n + 1) )
    ans <- c(m1, m2)
    return( ans )
  }
  lM2 <- lapply( 1:length(lE), FUN = FUN )
  m1 <- unlist( lapply( lM2, FUN = function(x) { x[1] } ) )
  m2 <- unlist( lapply( lM2, FUN = function(x) { x[2] } ) )
  
  boxplot(m1)
  boxplot(m2)
    
  
  
  
  
  
  
  
  