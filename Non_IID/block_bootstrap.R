  
  
  ############################################################################
  ### BLOCK BOOTSTRAP METHODS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.04.2021
  # First version:    09.04.2021
  # --------------------------------------------------------------------------
  
  
  # Notes:
  # How does the grid of over the probability simplex look like in the 
  # circular block bootstrap? First, the dimension of the simplex increases
  # to account for the circularity. As with the ordinary bootstrap, each sampling
  # defines a weight vector (point on the simplex). However, with the block bootstrap,
  # elements of the weights vector have to be consecutive. E.g., the centroid is not 
  # allowed.
  # Probably, the block bootstrap grid nodes are not 'typical'.
  
  #// Idea: use block bootstrap to compute distribution of Mahalanobis distance (turbulence).

  #// Idea: use sinus wave as conditioning variable, compute conditional density
  #         for a grid of conditioning values. Can this replace the block bootstrap?
  
  
  require(DAARC)
  
  Y <- log(1 + GPO::Data[ ,1])
  z <- scale(Y, TRUE, FALSE)^2
  n_sim <- 10^3
  block_length <- 20
  
  
  
  # --------------------------------------------------------------------------
  # Using package boot
  # --------------------------------------------------------------------------
 
  # ends <- cbind(i.a$starts[r, ], i.a$lengths)
  # inds <- apply(ends, 1L, make.ends, n)
  # inds <- if (is.list(inds)) 
  #   matrix(unlist(inds)[1L:n.sim], n.sim, 1L)
  # else matrix(inds, n.sim, 1L)
  
  uniques <- unique(inds)
  tmp <- unlist( lapply( uniques, FUN = function(x) { sum( inds %in% x ) } ) )
  w <- inds * 0 
  w[uniques, ] <- tmp
  w <- w / sum(w)
  plot(w)
  
    
  # Ordinary bootstrap
  B <- boot( data = z, 
             statistic = meanBoot,
             R = n_sim,
             sim = "ordinary" )
  
  # Fixed block bootstrap with length block_length
  # debugonce( tsboot )
  # debuognce( boot:::ts.array )
  FB <- tsboot( tseries = z, 
                statistic = mean,
                R = n_sim, 
                l = block_length, 
                sim = "fixed" )
  
  # Stationary bootstrap with mean block length block_length
  SB <- tsboot( tseries = z, 
                statistic = mean,
                R = n_sim, 
                l = block_length, 
                sim = "geom" )
                
 
  theta_mat <- cbind( ordinary = B$t,
                      block = FB$t,
                      stationary = SB$t )
  ldens <- apply( theta_mat, 2, density )
  plot.ldensity( ldens )
  abline( v = B$t0 )
  abline( v = FB$t0, col = 2 )
  abline( v = SB$t0, col = 3 )
  
  
  
  edit( boot:::ts.array )
  
  edit( tsboot )
  edit( boot.array )
  
  A <- boot.array( FB )
  str(A)
  headleft(A)
  
  
  ###
  Idx <- boot:::ts.array( n = nrow(Y),
                          n.sim = nrow(Y),
                          R = n_sim,
                          l = block_length,
                          sim = "fixed", 
                          endcorr = FALSE )
  
  
  headleft(Idx$starts)
  dim(Idx$starts)
  
  
  
  
  # --------------------------------------------------------------------------
  # Using package tseries
  # Note: sampling is done in C (see tseries:::boot.sample)
  # --------------------------------------------------------------------------
  
  require(tseries)
  
  # Generate AR(1) process
  n <- 500  
  a <- 0.6
  e <- rnorm(n+100)  
  x <- double(n+100)
  x[1] <- rnorm(1)
  for(i in 2:(n+100)) {
    x[i] <- a * x[i-1] + e[i]
  }
  x <- ts(x[-(1:100)])
  
  
  # debugonce( tsbootstrap )
  BB <- tsbootstrap( x = x,
                     nb = n_sim,
                     type = "block",
                     statistic = mean,
                     m = 1,
                     b = block_length )
  str(BB)

  
  SB <- tsbootstrap( x = Y,
                     nb = n_sim,
                     statistic = meanBoot,
                     m = 1,
                     b = block_length )
  
  edit( tsbootstrap )
  
  
  
  
  # --------------------------------------------------------------------------
  # Instead of block bootstrap, do ordinary bootstrap on rolling means
  # --------------------------------------------------------------------------
  
  z_roll <- applyRoll( Data = z,
                       Width = block_length,
                       By = 1,
                       FUN = mean )
  
  plot( cbind(z, z_roll), plot.type = "single" )

  
  meanBootCustom <- function(x, idx)
  {
    m <- 45
    idx <- sample(idx)[1:m]
    if ( NCOL(x) > 1 ) {
      ans <- apply(x[idx, ], 2, mean)
    } else {
      ans <- mean(x[idx])
    }
    return( ans )
  }
  
  meanBootCustomBayes <- function(x, idx)
  {
    m <- 45
    idx <- sample(idx)[1:m]
    w <- rdirichlet(n = 1, rep(1, m))
    if ( NCOL(x) > 1 ) {
      ans <- apply(x[idx, ], 2, function(y) { sum(w * y) } )
    } else {
      ans <- sum(w * x[idx])
    }
    return( ans )
  }

  # Ordinary bootstrap
  B1 <- boot( data = z, 
             statistic = meanBoot,
             # statistic = meanBootCustom,
             R = n_sim,
             sim = "ordinary" )
  
  B2 <- boot( data = z_roll, 
              statistic = meanBootCustom,
              R = n_sim,
              sim = "ordinary" )
  
  # debugonce( meanBootCustomBayes )
  B3 <- boot( data = z_roll, 
              statistic = meanBootCustomBayes,
              R = n_sim,
              sim = "ordinary" )
  
  # Fixed block bootstrap
  FB <- tsboot( tseries = z, 
                statistic = mean,
                R = n_sim, 
                l = block_length, 
                sim = "fixed" )
  
  # Stationary block bootstrap
  SB <- tsboot( tseries = z, 
                statistic = mean,
                R = n_sim, 
                l = block_length, 
                sim = "geom" )
  
  
  # Plot
  theta_mat <- cbind( ordinary = as.numeric(B1$t),
                      o2 = as.numeric(B2$t),
                      o2_bayes = as.numeric(B3$t),
                      block = as.numeric(FB$t),
                      stationary = as.numeric(SB$t) )
  ldens <- apply( theta_mat, 2, density )
  plot.ldensity( ldens, colors = 1:ncol(theta_mat), fillin = FALSE, main = "Bootstrap distributions" )
  abline( v = B1$t0 )
  abline( v = B2$t0, col = 2 )
  abline( v = B3$t0, col = 3 )
  abline( v = FB$t0, col = 4 )
  abline( v = SB$t0, col = 5 )
  
 
  mean(z); mean(z_roll)
  B1$t0; B2$t0; B3$t0; FB$t0; SB$t0
  
  
  
  
  
  
  
  
  
  
