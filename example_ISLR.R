  
  
  ############################################################################
  ### EXAMPLE ISLR p.187
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     26.07.2021
  # First version:    26.07.2021
  # --------------------------------------------------------------------------
  
  
  #// Example fails because of division (nonlinear statistic)
  
  
  
  require(GPO)
  require(boot)
  
  n <- 100
  s11 <- 1
  s22 <- 1.25
  s12 <- 0.5
  covmat <- matrix( c(s11, s12, s12, s22), nrow = 2, ncol = 2, byrow = TRUE )
  mu <- c(0, 0)
  # X <- rmvnorm( n = n, mean = mu, sigma = covmat )
  X <- miscolz::rmvnorm( nObs = n, nAssets = 2, mu = mu, sigma = covmat )
  
  a_num <- (s22 - s12) 
  a_denom <- (s11 + s22 - 2 * s12)
  a <- a_num / a_denom
  a  
  
  cov(X)
  covmat
  
  Y <- scale( X, TRUE, FALSE )
  Z_num <- Y[ ,2]^2 - Y[ ,1] * Y[ ,2]
  Z_denom <- Y[ ,1]^2 + Y[ ,2]^2 - 2 * Y[ ,1] * Y[ ,2]
  
  mean(Z_num) / (1 - 1/n)
  a_num 
  
  mean(Z_denom) / (1 - 1/n)
  a_denom
  
  mean(Z_num) / mean(Z_denom)
  sum(Z_num) / sum(Z_denom)
  a
  
  Z <- Z_num * 1 / sum(Z_denom)
  sum(Z)  
  a
  
  var(Z) / n
  sd(Z) / sqrt(n) * n
  
  
  
    
  bootFUN <- function(x, idx)
  {
    covmat <- cov(x[idx, ])
    ans <- (covmat[2, 2] - covmat[1, 2]) / (covmat[1, 1] + covmat[2, 2] - 2 * covmat[1, 2])
    return( ans )
  }
  Boot <- boot::boot( data = X, 
                      statistic = bootFUN,
                      R = 10^3 )
  Boot2 <- boot::boot( data = Z,
                       statistic = meanBoot,
                       R = 10^3 )
  
  mean(Boot$t)
  sd(Boot$t)
  var(Boot$t)
  
  sd(Boot$t); sd(Boot2$t)

  
  
  # Sharpe ratio
  x <- rnorm(100, 0.1, 2)
  y <- (x - mean(x))^2
  mean(x) / sd(x)  
  mean(x) / sqrt(mean(y))    
  mean( x / sqrt(mean(y)) )
    
  
  
  
  
  
  
  