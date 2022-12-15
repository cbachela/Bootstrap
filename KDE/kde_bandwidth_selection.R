  
  
  ############################################################################
  ### WEIGHTED KERNEL DENSITY ESTIMATION - BANDWIDTH SELECTION
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
  # # x <- runif(100)
  # x_eval <- seq(from = 0, to = 1, length.out = 10^2+1)
  x <- abs(rnorm(100))^rnorm(100)
  x_eval <- quantile( x, seq(from = 0, to = 1, length.out = 10^2+1) )
  n <- length(x)
  
  dens <- density(x, bw = "SJ")
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
  # Same example as above but using random weighting
  # --------------------------------------------------------------------------
  
  n_sim <- 10^3
  n <- length(x)
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  M <- (samples / h) %*% mat_h # * n
  
  plot( x = x_eval, y = M[1, ], ylim = range(M) )
  for ( i in 1:nrow(M) ) {
    lines( x = x_eval, y = M[i, ], col = 4 )
  }
  for ( i in 1:nrow(mat_h) ) {
    lines( x = x_eval, y = mat_h[i, ], col = 2 )
  }
  abline( v = x, col = "darkgrey" )
  lines( dens )
  
  m2 <- apply( M, 2, var )
  plot( m2, type = "o" )
  
  
  
  f_h <- apply(mat_h, 2, sum)
  tmp <- apply( M, 1, function(f) { (f - f_h)^2 } )
  mise <- mean( tmp )
  mise
  
  
  # -------------------
  # Find h that minimizes mise
  # -------------------
  
  lmat <- list()
  h_vec <- seq( from = h/10, to = h*10, length.out = 100 )
  for ( k in seq(along = h_vec) ) {
    h_tmp <- h_vec[k]
    mat <- matrix( 0, nrow = length(x), ncol = length(x_eval) )
    for ( i in 1:nrow(mat) ) {
      for ( j in 1:ncol(mat) ) {
        mat[i, j] <- kernelFunction( x = (x_eval[j] - x[i]) / h_tmp, 
                                     method = "gaussian" )
      }
    }
    lmat[[k]] <- mat
  }
  
  FUN <- function(i)
  {
    M <- samples %*% (lmat[[i]] / h_vec[i])
    tmp <- apply( M, 1, function(f) { mean((f - f_h)^2) } )
    mise <- sum( tmp )
    # m2 <- apply( M, 2, var )
    # mise = sum(m2)
    return( mise )
  }
  # f_h <- apply( mat_h, 2, sum ) / (h*n)
  k <- 100; f_h <- apply( lmat[[k]], 2, sum ) / (h_vec[k]*n)
  mise <- unlist( lapply( 1:length(lmat), FUN = FUN ) )
  
  plot( x = h_vec, y = mise, type = "o" )
  
  h
  h_vec[ which(mise == min(mise)) ]
  
  
  density(x)$bw
  density(x, bw = "SJ")$bw
  
  
  
  # -------------------
  # Classical bootstrap
  # -------------------
  
  x_boot <- matrix( 0, nrow = n_sim, ncol = n )
  for ( i in 1:n_sim ) {
    idx <- sample( x = 1:n, size = n, replace = TRUE )
    x_boot[i, ] <- x[idx]
  }
  
  
  lmat <- list()
  h_vec <- seq( from = h/10, to = h*10, length.out = 100 )
  for ( k in seq(along = h_vec) ) {
    h_tmp <- h_vec[k]
    mat <- array( 0, dim = c(nrow(x_boot), nrow = length(x), ncol = length(x_eval)) )
    for ( i in 1:dim(mat)[[2]] ) {
      for ( j in 1:dim(mat)[[3]] ) {
        for ( s in 1:dim(mat)[[1]] ) {
          mat[s, i, j] <- kernelFunction( x = (x_eval[j] - x_boot[s, i]) / h_tmp, 
                                          method = "gaussian" )
        }
      }
    }
    lmat[[k]] <- mat
  }
  
  
  FUN <- function(i)
  {
    fmat <- matrix( NA, nrow = dim(lmat[[i]])[1], ncol = dim(lmat[[i]])[3] )
    for( s in 1:nrow(fmat) ) {
      fmat[s, ] <- apply( lmat[[i]][s,, ], 2, sum ) / (h_vec[i]*n)
    }
    # tmp <- scale( fmat, TRUE, FALSE )^2
    tmp <- t( apply( fmat, 1, function(f) { (f - f_h)^2 } ) )
    mise <- mean( apply( tmp, 1, sum) )
    return( mise )
  }
  
  # f_h <- apply( mat_h, 2, sum ) / (h*n)
  f_h <- apply( lmat[[1]], 2, sum ) / (h*n)
  mise <- unlist( lapply( 1:length(lmat), FUN = FUN ) )
  
  plot( x = h_vec, y = mise, type = "o" )
  
  h
  h_vec[ which(mise == min(mise)) ]
  
  
  
  
  
  
