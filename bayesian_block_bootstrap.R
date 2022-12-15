  
  
  ############################################################################
  ### BAYESIAN BLOCK BOOTSTRAP METHODS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.08.2021
  # First version:    16.08.2021
  # --------------------------------------------------------------------------
  
  
  # We try to find a Bayesian analog to the ordinary block bootstrap and it's 
  # exact solution. I.e., we need an analog to the discrete m out of n, m < n sampling 
  # for the continuous Dirichlet sampling.
  
  # Question:
  # Is the expected 2-norm the same for n dim with zeros in the Dirichlet parameter vs m dim?
  
  
  
  require(DAARC)
  
  B <- 10^4 + 1
  n_boot <- 50
  Y <- log(1 + GPO::Data[ ,1])
  z <- scale(Y, TRUE, FALSE)^2
  n_sim <- 10^3
  n <- nrow(z)
  block_length <- 20
  
  
  
  
  
  z_roll <- applyRoll( Data = z,
                       Width = block_length,
                       By = 1,
                       FUN = mean )
  
  plot( cbind(z, z_roll), plot.type = "single" )
  
  
  
  
  m <- 45 # =~ n / block_length
  alpha_n <- rep(1/n, n) * n
  alpha_m <- c( rep(1 / m, m), rep(0, n-m) ) * m   # ~~~~~~~~~~~~~
  sum(alpha_n); sum(alpha_m)
  
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha_n) - a, m = 2 ) }
  beta_moments_n <- unlist( lapply( alpha_n, FUN = betaFUN ) )
  beta_moments_m <- unlist( lapply( alpha_m, FUN = betaFUN ) )
  
  sum(beta_moments_n)
  sum(beta_moments_m) # different
  sum(beta_moments_m) * n/m # same same
  
  
  cbind( beta_moments_n, 
         beta_moments_m, 
         (alpha_n^2 + alpha_n) / (n^2 + n) )
  
  
  
  bb_mat_n <- matrix( NA, nrow = B, ncol = n_boot )
  bb_mat_m <- matrix( NA, nrow = B, ncol = n_boot )
  bb_mat_mn <- matrix( NA, nrow = B, ncol = n_boot )
  lP_n <- lP_m <- lP_mn <- list()
  
  for ( j in 1:n_boot ) {
    
    lP_n[[j]] <- rdirichlet( n = B, alpha = alpha_n )
    bb_mat_n[ ,j] <- apply( lP_n[[j]], 1, function(p) { sum(p * z) } )
    
    lP_m[[j]] <- rdirichlet( n = B, alpha = rep(1, m) )
    z_m <- sample(z)[1:m]
    bb_mat_m[ ,j] <- apply( lP_m[[j]], 1, function(p) { sum(p * z_m) } )
    
    lP_mn[[j]] <- rdirichlet( n = B, alpha = sample(alpha_m) )
    bb_mat_mn[ ,j] <- apply( lP_mn[[j]], 1, function(p) { sum(p * z) } )
    
  }
  
  
  do.call( cbind, lapply( lP_n, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) )
  do.call( cbind, lapply( lP_m, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) )
  do.call( cbind, lapply( lP_mn, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) )
  
  range( do.call( cbind, lapply( lP_n, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) ) )
  range( do.call( cbind, lapply( lP_m, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) ) )
  range( do.call( cbind, lapply( lP_mn, FUN = function(P) { apply(P, 1, lpNorm, p = 2 )} ) ) )
  
  
  
  ldens_bb <- list( n = density(as.numeric(bb_mat_n)),
                    m = density(as.numeric(bb_mat_m)),
                    mn = density(as.numeric(bb_mat_mn)) )
  plot.ldensity( ldens_bb, colors = 1:length(ldens_bb), fillin = FALSE )
  
  
  
  
  # Compare to results from script block_bootstrap.R
  ldens <- c( ldens_bb, ldens_boot )
  plot.ldensity( ldens, colors = 1:length(ldens), fillin = FALSE )
  
  
  
  
  
  
  

  
  
  
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
    idx <- sample(idx)[1:45]
    # idx <- sample(idx)[1:100]
    # idx <- sample(idx)[1:35]
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
                      # o2_bayes = as.numeric(B3$t),
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
  
  
  
  
  
  
  
  
  
  
