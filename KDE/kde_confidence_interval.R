  
  
  ############################################################################
  ### WEIGHTED KERNEL DENSITY ESTIMATION - CONFIDENCE INTERVAL
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     27.05.2021
  # First version:    27.05.2020
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
     
  
  
  # --------------------------------------------------------------------------
  # 1D illustration of how KDE is constructed
  # --------------------------------------------------------------------------
  
  # x <- c(0.18, 0.2, 0.27, 0.6, 0.8)
  # x_eval <- quantile( x, seq(from = 0, to = 1, length.out = 10^2+1) )
  # x <- runif(100)
  # x_eval <- seq(from = 0, to = 1, length.out = 10^2+1)
  # x <- abs(rnorm(100))^rnorm(100)
  # x_eval <- quantile( x, seq(from = 0, to = 1, length.out = 10^2+1) )
  x <- seq(from = 1010, to = 1070, by = 10)
  x <- c(x, rep(1, 81 - length(x)))
  x_eval <- quantile( x, seq(from = 0, to = 1, length.out = 10^2+1) )
  
  n <- length(x)
  # dens <- density(x, bw = "SJ")
  dens <- density(x)
  h <- dens$bw
  mat_h <- matrix( 0, nrow = length(x), ncol = length(x_eval) )
  for ( i in 1:nrow(mat_h) ) {
    for ( j in 1:ncol(mat_h) ) {
      mat_h[i, j] <- kernelFunction( x = (x_eval[j] - x[i]) / h, 
                                     method = "gaussian" )
    }
  }  
  
  # mat_h <- t( apply( mat_h, 1, function(x) { x / sum(x) } ) )  # make it right-stochastic
  # apply( mat_h, 1, sum )
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), col = 4, type = "l" )
  abline( v = x, col = "darkgrey" )
  for ( i in 1:nrow(mat_h) ) {
    lines( x = x_eval, y = mat_h[i, ], col = 2 )
  }
  lines(dens)
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Confidence interval based on classical bootstrap
  # --------------------------------------------------------------------------
  
  n_sim <- 10^3
  dens <- density(x, bw = "SJ")
  h <- dens$bw
  
  fmat_boot <- matrix( NA, nrow = n_sim, ncol  = length(x_eval) )
  for ( s in 1:n_sim ) {
    
    idx <- sample( x = 1:n, size = n, replace = TRUE )
    x_boot <- x[idx]
    
    mat <- matrix( 0, nrow = length(x), ncol = length(x_eval) )
    for ( i in 1:nrow(mat_h) ) {
      for ( j in 1:ncol(mat_h) ) {
        mat[i, j] <- kernelFunction( x = (x_eval[j] - x_boot[i]) / h, 
                                     method = "gaussian" )
      }
    } 
    fmat_boot[s, ] <- apply( mat, 2, sum ) / (h*n)
  }
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), 
        ylim = range(fmat), col = 2, type = "l" )
  abline( v = x, col = "darkgrey" )
  for ( s in 1:nrow(fmat_boot) ) {
    lines( x = x_eval, y = fmat_boot[s, ], col = 4 )
  }
  for ( i in 1:nrow(mat_h) ) {
    lines( x = x_eval, y = mat_h[i, ], col = 2 )
  }
  lines(dens)
  
  
  
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), 
        ylim = range(fmat), col = 2, type = "l" )
  lines( x = x_eval, y = apply( fmat_boot, 2, quantile, 0.025), col = 4 )
  lines( x = x_eval, y = apply( fmat_boot, 2, quantile, 0.975), col = 4 )
  
  
  
  
  # --------------------------------------------------------------------------
  # Confidence interval based on weighted bootstrap
  # --------------------------------------------------------------------------
  
  n_sim <- 10^3
  dens <- density(x, bw = "SJ")Q
  h <- dens$bw
  
  alpha <- rep(1, n) 
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  fmat <- (samples %*% mat_h / h) 
  
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), 
        ylim = range(fmat), col = 2, type = "l" )
  abline( v = x, col = "darkgrey" )
  for ( s in 1:nrow(fmat) ) {
    lines( x = x_eval, y = fmat[s, ], col = 4 )
  }
  for ( i in 1:nrow(mat_h) ) {
    lines( x = x_eval, y = mat_h[i, ], col = 2 )
  }
  lines(dens)
  
  
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), 
        ylim = range(fmat), col = 2, type = "l" )
  lines( x = x_eval, y = apply( fmat, 2, quantile, 0.025), col = 4 )
  lines( x = x_eval, y = apply( fmat, 2, quantile, 0.975), col = 4 )
  
  
  
  # --------------------------------------------------------------------------
  # Confidence interval using Varsi
  # --------------------------------------------------------------------------
  
  FUN <- function( x, th = 0.025 )
  {
    varsiFindQuantile( z = x, th = th )
  }
  
  lb <- apply( mat_h / h, 2, FUN = FUN, th = 0.025 )
  ub <- apply( mat_h / h, 2, FUN = FUN, th = 0.975 )
  
  
  plot( x = x_eval, y = apply( mat_h, 2, sum ) / (h*n), 
        ylim = range(fmat), col = 1, type = "l" )
  lines( x = x_eval, y = apply( fmat, 2, quantile, 0.025), col = 2 )
  lines( x = x_eval, y = apply( fmat, 2, quantile, 0.975), col = 2 )
  lines( x = x_eval, y = apply( fmat_boot, 2, quantile, 0.025), col = 3 )
  lines( x = x_eval, y = apply( fmat_boot, 2, quantile, 0.975), col = 3 )
  lines( x = x_eval, y = lb, col = 4 )
  lines( x = x_eval, y = ub, col = 4 )
  
  
  
  
  
  
  