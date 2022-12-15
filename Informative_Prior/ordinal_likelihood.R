  
  
  ############################################################################
  ### ORDINAL LIKELIHOOD
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.05.2021
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  # Ordering view in the parameter of the multinomial distribution:
  # p_1 >= p_2 >= ... >= p_n
  # I.e., orderging of the \alpha parameter of the Dirichlet distribution
  
  # Research question:
  # Can we bound the variation of the statistic like in the exact boot se?
  # Reduces to the question whether we can bound 2-norm of the samples over the 
  # constrained simplex.
  
  # Compare uniform samples over constrained simplex to 
  # Dirichlet samples with \alpha = E(P)
  
  
  
  # Sections:
  # 
  # Specifications
  # Geometric random walk sampling
  # Exponentially weighted average
  # Shadow Dirichlet
  # Acceptance-rejection of Gridnodes with probability rescaling
  
  
  
  
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
  source( "H:/R/RP/R/class_Polytope.R" )
  
  
  
  # --------------------------------------------------------------------------
  # Specifications
  # --------------------------------------------------------------------------
  
  
  n <- 3
  n_sim <- 10^5
  set.seed(1111)
  x <- rnorm(n)
  z <- (x - mean(x))^2
  # names(z) <- aviationAlphabet()[1:n]
  names(z) <- paste0("X", 1:n)
  
  
  # --------------------------------------------------------------------------
  # Geometric random walk sampling
  # --------------------------------------------------------------------------
  
  Constraints <- constraints( selection = names(z) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint( lower = setNames(rep(0, n), names(z)),
                                               upper = setNames(rep(1, n), names(z)) )
  addConstraint(Constraints) <- orderingConstraint( Names = names(z),
                                                    ordering = 1:n,
                                                    eps = 1e-04 )
  P_tmp <- as.polytope(Constraints) 

  P <- Polytope$new( A = P_tmp@A,
                     b = P_tmp@b, 
                     sense = P_tmp@sense )

  interior_point <- setNames( (1:n) / sum(1:n), names(z) )
  PS <- PolytopeSampler$new( polytope = P )
  PS$setCtrl( algo = "hitnrun",
              interior_point = interior_point )
  # PS$polytope$isInterior( x = interior_point )
  PS$transform()
  lSamples <- PS$polytope$sample( algo = "hitnrun",
                                  n_sim = n_sim,
                                  x0 = PS$spec$interior_point )
  PS$samples <- lSamples
  PS$transformBack()
  PS$samples
  samples <- PS$samples$S
  
  
  lp_samples <- apply( samples, 1, lpNorm, p = 2 )
  plot( density( lp_samples^2 ) )

  
  # Use Dirichlet
  alpha <- apply(samples, 2, mean)
  Omega <- rdirichlet( n = n_sim,
                       alpha = alpha * n )
  lp_omega <- apply( Omega, 1, lpNorm, p = 2 )
  
  cbind( alpha, apply(Omega, 2, mean) )

  
  plot( density(lp_samples) )  
  lines( density(lp_omega), col = 2 )
  abline( v = lpNorm(alpha, p = 2) )
  abline( v = lpNorm(apply(Omega, 2, mean), p = 2), col = 2 )
    
    
  
  
  
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  samples_simplex <- getSamples( rp( simplex( d = n ), rpCtrl( n_sim = 10^4 ) ) )
  
  plot3d( x = samples_simplex[ ,1], y = samples_simplex[ ,2], z = samples_simplex[ ,3], 
          xlab = names(z)[1], ylab = names(z)[2], zlab = names(z)[3], col = "lightgrey" )
  points3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "darkgrey" )
  points3d( x = Omega[ ,1], y = Omega[ ,2], z = Omega[ ,3], col = "green" )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Exponentially weighted average
  # --------------------------------------------------------------------------
  
  B <- 10^4
  tau <- 52
  X <- log( 1 + GPO::Data[ ,1] )
  X <- head(X, 300)
  mu <- mean(X)
  mu_ewma <- mean( weightEWMA( X = X, tau = tau ) )
  n <- nrow(X)
  z <- setNames( as.numeric(X), rownames(X) ) 
  
  lambda <- exp(-log(2) / tau)
  wt <- lambda^(0:(n-1))
  wt <- rev( length(wt) * wt / sum(wt) ) / n
  sum( wt * z )
  
  
  tau <- 10
  lambda <- exp(-log(2) / tau)
  wt <- lambda^(0:(n-1))
  wt <- rev( length(wt) * wt / sum(wt) ) / n
  plot( wt )
  mu_ewma <- sum( wt * z )
  
  
  # Weighted ordinary bootstrap
  Boot <- boot::boot( data = z,
                      statistic = meanBoot,
                      R = B,
                      weights = wt )
  
 
  
  # Using Volesti
  Constraints <- constraints( selection = names(z) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint( lower = setNames(rep(0, n), names(z)),
                                               upper = setNames(rep(1, n), names(z)) )
  addConstraint(Constraints) <- orderingConstraint( Names = names(z),
                                                    ordering = 1:n, 
                                                    eps = 0 )
  P <- as.polytope(Constraints) 
  
  # Interior point used for starting the MCMC
  x0 <- getInteriorPoint(P)
  x0 <- wt
  
  # isInterior(P, x0)
  # isInterior(P, wt)
  # head( cbind( P@A %*% x0, P@b ) )
  # tail( cbind( P@A %*% x0, P@b ) )
  # idx <- ( P@A %*% x0 <= P@b )
  # which(idx == 0)
  # headleft(P@A[which(idx == 0), ])
  # head(P@sense[which(idx == 0)])
  # head( cbind( P@A %*% x0, P@b )[which(idx == FALSE), ] )
  
  # Specifications
  rp_ctrl <- rpCtrl( n_sim = 10^4,
                     thin = 10^3,
                     # algo = "billiard",
                     algo = "volesti",
                     volesti = list( density = "uniform",
                                     walk = "john" ),
                     x0 = x0 )
                     # x0 = NULL )
  
  # debugonce( .rpP )
  tic <- Sys.time()
  RPP <- rp( object = P, spec = rp_ctrl )
  (toc <- Sys.time() - tic)
  samples <- getSamples(RPP)
  theta_volesti <- apply( samples, 1, function(p) { sum( p * z ) } )
  lp_samples <- apply( samples, 1, lpNorm, p = 2 )
  
  
  
  # Specifications
  rp_ctrl <- rpCtrl( n_sim = 10^4,
                     thin = 10^3,
                     x0 = x0 )
  # debugonce( .rpP )
  tic <- Sys.time()
  RPP <- rp( object = P, spec = rp_ctrl )
  (toc <- Sys.time() - tic)
  
  
  
  
  plot( samples[ ,ncol(samples)] )
  
  
  
  
  
  
  
  # Use Dirichlet
  alpha <- apply(samples, 2, mean)
  Omega <- rdirichlet( n = rp_ctrl$n_sim,
                       alpha = alpha * n )
  theta_bb <- apply( Omega, 1, function(p) { sum(p * z) } )
  
  lp_omega <- apply( Omega, 1, lpNorm, p = 2 )
  
  
  
  plot( x = wt, y = apply(samples, 2, mean) )
  plot( x = wt, y = apply(Omega, 2, mean) )
  
  
  plot( wt )
  lines( samples[1, ] )
  lines( Omega[1, ] )
  
  
  ldens <- list( boot = density(Boot$t),
                 volesti = density(theta_volesti),
                 bb = density(theta_bb) )
  plot.ldensity( ldens, fillin = FALSE, col = 1:3 )
  abline( v = mu )
  abline( v = mu_ewma, col = 2 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Shadow Dirichlet
  # Bela A. Frigyik, Maya R. Gupta, and Yihua Chen
  # --------------------------------------------------------------------------
  
  # Findings:
  
  # - Cannot just use betaMoments on transformed alpha bcs shadow Dirchlet
  #   density has an additional scaling term det(|M|).
  # - Identity $SE(\hat{\theta}) = s^2 ( \mathds{E}(||\omega||_2^2) - ||\bar{\omega}||_2^2 )$ 
  #   does not hold when the domain is not the entire simplex as in the shadow case.
  # - Uniform samples on the simplex and on the ordinal constrained simplex have the 
  #   same expected squared 2-norm. 
  # - Can we use the expected 2-norm towards the center of mass?
  # - What if we shift the new domain to the center of the simplex?
  # - Does that mean that samples on the constrained simplex to not concentrate 
  #   near the boundary? 
  # - Would we get the same result when sampling the consrained simplex with a 
  #   random walk rather than doing the shadow Dirichlet transformation?
  
  #// For p = 2,
  #   we observe another peculiar phenomena that each vector's distance to the centroid is
  #   approximately equal to the l2-norm of the centroid, i.e., ||xn k ??? x¯n||2 ??? ||x¯n||2 ??? ??pn.
  #   --> does that mean \mathds{E}(||\omega||_2^2) - ||\omega||_2^2 = \mathds{E}(||\bar{\omega}||_2^2) ?
  
  
  
  
  require(rgl)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  source( "H:/R/RP/R/class_Polytope.R" )
  source( "H:/R/RP/R/class_Simplex.R" )
  source( "H:/R/RP/R/class_Ellipsoid.R" )
  source( "H:/R/RP/R/class_Sampler.R" )

  
  
  
  
  
  n <- 3
  n_sim <- 10^5
  set.seed(9999)
  x <- rnorm(n)
  z <- (x - mean(x))^2 * 10^3
  # z <- rev(sort(z))
  # names(z) <- aviationAlphabet()[1:n]
  names(z) <- paste0("X", 1:n)
  w_bar <- rep(1/n, n)
  # w_bar <- (1:n / sum(1:n))
  
  
  S <- Simplex$new( d = n )
  S
  samples <- S$runif( n_sim = n_sim )
  
  # Or use Dirichlet
  alpha <- w_bar * n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  
  
  # Generate the transformation matrix
  
  # tic <- Sys.time()
  # M <- diag(n) * 0
  # for ( k in 1:ncol(M) ) {
  #   M[ ,k] <- c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) )
  # }
  # (toc1 <- Sys.time() - tic)
  
  tic <- Sys.time()
  M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
  (toc2 <- Sys.time() - tic)
  
  # ###
  # Sigma <- cov(samples)
  # M <- t( eigen(Sigma)$vectors )
  # ###
  
  
  
  
  
  # Transform samples
  shadow <- t( apply( samples, 1, function(p) { M %*% matrix(p, ncol = 1) } ) )
  w_shadow_bar <- M %*% matrix(w_bar, ncol = 1)
  cbind( w_shadow_bar, apply( shadow, 2, mean) )
  
  # Shift transformed samples to have the origin as centroid
  shadow_shifted <- t( apply( shadow, 1, function(x) { x + as.numeric(w_bar - w_shadow_bar) } ) )
  apply( shadow_shifted, 2, mean )
  

  # Transform Dirichlet parameter
  alpha_prime <- M %*% alpha
  cbind( alpha_prime / n, M %*% alpha/n, w_shadow_bar ) # same same
  
  
  # FEV
  samples_fev <- na.omit( fevBias( x = samples, q = 1000 ) )
  # Transformed fev-biased samples
  shadow_fev <- t( M %*% t(samples_fev) )
  
  
  # Sample from Dirichlet with transformed parameter
  # Note: contrary to shadow Dirichlet, the domain is still the entire simplex.
  P <- rdirichlet( n = nrow(samples), alpha = alpha_prime )
  
  # Sample from sphere with radius E(lp_samples)
  E <- Ellipsoid$new( centre = rep(0, n),
                      shape = diag(n),
                      rhs = mean( apply( samples, 1, lpNorm, p = 2 )^2 ) )
  ES <- EllipsoidSampler$new( ellipsoid = E )
  ES$setCtrl( n_sim = n_sim * 10,
              algo = "sphere",
              b_pushy = FALSE )
  ES$sample()
  samples_E <- ES$samples$S
  
  
  
  
  # Plot
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1 )
  points3d( x = P[ ,1], y = P[ ,2], z = P[ ,3], col = "green", size = 2 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 3 )
  points3d( x = z_shadow_bar[1], y = z_shadow_bar[2], z = z_shadow_bar[3], size = 10 )
  # points3d( x = samples_E[ ,1], y = samples_E[ ,2], z = samples_E[ ,3], col = 1, size = 5 )
  # points3d( x = shadow_shifted[ ,1], y = shadow_shifted[ ,2], z = shadow_shifted[ ,3], col = "red", size = 3 )
  points3d( x = samples_fev[ ,1], y = samples_fev[ ,2], z = samples_fev[ ,3], col = "red", size = 3 )
  points3d( x = shadow_fev[ ,1], y = shadow_fev[ ,2], z = shadow_fev[ ,3], col = "red", size = 5 )
  
  
  
  # Checks
  cbind( apply( samples, 2, mean ) * n, alpha )
  cbind( apply( shadow, 2, mean ) * n, alpha_prime, apply( P, 2, mean ) * n )
  
  
  
  # --------------------------------------------------------------------------
  # Next, compute expected 2-norm of shadow dirichlet pmf's 
  # and check whether standard error can still be computed exactly.
  # --------------------------------------------------------------------------
  
  
  lp_samples <- apply( samples, 1, lpNorm, p = 2 )
  lp_samples_rel <- apply( samples, 1, lpNorm, p = 2, x_bar = w_bar )
  lp_shadow <- apply( shadow, 1, lpNorm, p = 2 )
  lp_shadow_rel <- apply( shadow, 1, lpNorm, p = 2, x_bar = w_shadow_bar )
  lp_P <- apply( P, 1, lpNorm, p = 2 )
  lp_shadow_shifted <- apply( shadow_shifted, 1, lpNorm, p = 2 )
  
  
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  delta <- abs( lp_samples^2 - sum(beta_moments) )
  idx <- which(delta < 1e-03)
  
  # Plot
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 2 )
  points3d( x = samples[idx, 1], y = samples[idx, 2], z = samples[idx, 3], col = "darkgrey", size = 3 )
  points3d( x = shadow[idx, 1], y = shadow[idx, 2], z = shadow[idx, 3], col = "red", size = 3 )
  points3d( x = w_bar[1], y = w_bar[2], z = w_bar[3], size = 5, col = "darkgrey" )
  points3d( x = w_shadow_bar[1], y = w_shadow_bar[2], z = w_shadow_bar[3], size = 5, col = "red" )
  
  
  
  
  
  plot( lp_samples, lp_shadow)
  plot( lp_samples_rel, lp_shadow_rel)
  
  mean( lp_samples^2 ); mean( lp_shadow^2 ); mean( lp_shadow_shifted^2 ); mean( lp_P^2 )
  mean( lp_samples_rel^2 ); mean( lp_shadow_rel^2 )
  boxplot( cbind(lp_samples, lp_shadow, lp_P) )
  
  
  
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  betaFUNPrime <- function(a) { betaMoments( a = a, b = sum(alpha_prime) - a, m = 2 ) }
  beta_moments_prime <- unlist( lapply( alpha_prime, FUN = betaFUNPrime ) )
  
  
  sum(beta_moments); mean(lp_samples^2)
  sum(beta_moments_prime); mean(lp_shadow^2); mean(lp_P^2)
  
  cbind( beta_moments, apply( samples, 2, function(x) { mean(x^2) } ) )
  cbind( beta_moments_prime, apply( P, 2, function(x) { mean(x^2) } ) )
  cbind( beta_moments_prime, apply( shadow, 2, function(x) { mean(x^2) } ) )
  
  
  
  
  theta <- apply( samples, 1, function(p) { sum( p * z ) } )
  theta_shadow <- apply( shadow, 1, function(p) { sum( p * z ) } )
  theta_shadow_shifted <- apply( shadow_shifted, 1, function(p) { sum( p * z ) } )
  theta_P <- apply( P, 1, function(p) { sum( p * z ) } )
  theta_fev <- apply( samples_fev, 1, function(p) { sum(p * z) } )
  theta_shadow_fev <- apply( shadow_fev, 1, function(p) { sum(p * z) } )
  
  
  boxplot( cbind(theta, theta_shadow, theta_P, theta_shadow_shifted) )
  
  
  par(mfrow = c(2, 2))
  plot( density(theta) ) 
  plot( density(theta_fev) ) 
  abline( v = z )
  plot( density(theta_shadow) ) 
  plot( density(theta_shadow_fev) )
 
  plot( theta_shadow_fev )
  
  z
  var(z)
  var(theta_fev)
  t( M %*% diag(rep(1, n)) )
  apply( t( M %*% diag(rep(1, n)) ), 1, function(p) { sum(p * z) } )
  var( apply( t( M %*% diag(rep(1, n)) ), 1, function(p) { sum(p * z) } ) )
  
  
 
  
  alpha_0 <- sum(alpha)
  Ewi2 <- (alpha * (alpha + 1)) / (alpha_0 * (alpha_0 + 1))
  
  
  Ewij_emp <- matrix(n^n, nrow = n, ncol = n)
  Ewij_shadow_emp <- matrix(n^n, nrow = n, ncol = n)
  Ewij_P_emp <- matrix(n^n, nrow = n, ncol = n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      Ewij_emp[i, j] <- mean( samples[ ,i] * samples[ ,j] )
      Ewij_shadow_emp[i, j] <- mean( shadow[ ,i] * shadow[ ,j] )
      Ewij_P_emp[i, j] <- mean( P[ ,i] * P[ ,j] )
    }
  }
  
  Ewij_emp
  Ewij_shadow_emp
  Ewij_P_emp
  
  # Compare to analytical solution
  # ...
  
  Cov <- matrix(n^n, nrow = n, ncol = n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i == j ) {
        Cov[i, j] <- (alpha[i] * (alpha_0 - alpha[i])) / (alpha_0^2 * (alpha_0 + 1))
      } else {
        Cov[i, j] <- -(alpha[i] * alpha[j]) / (alpha_0^2 * (alpha_0 + 1))
      }
    }
  }
  
  
  alpha_0 <- sum(alpha)
  alpha_prime_0 <- sum(alpha_prime)
  Ewij <- matrix(n^n, nrow = n, ncol = n)
  Ewij_P <- matrix(n^n, nrow = n, ncol = n)
  Ewij_shadow <- matrix(n^n, nrow = n, ncol = n)
  for ( i in 1:n ) {
    for ( j in 1:n ) {
      if ( i == j ) {
        Ewij[i, j] <- (alpha[i] * (alpha[i] + 1)) / (alpha_0 * (alpha_0 + 1))
        Ewij_P[i, j] <- (alpha_prime[i] * (alpha_prime[i] + 1)) / (alpha_prime_0 * (alpha_prime_0 + 1))
        Ewij_shadow[i, j] <- (M %*% Cov %*% t(M))[i, j] + prod((M %*% alpha / alpha_0)[c(i, j)]) 
      } else {
        Ewij[i, j] <- (alpha[i] * alpha[j]) / (alpha_0 * (alpha_0 + 1))
        Ewij_P[i, j] <- (alpha_prime[i] * alpha_prime[j]) / (alpha_prime_0 * (alpha_prime_0 + 1))
        Ewij_shadow[i, j] <- (M %*% Cov %*% t(M))[i, j] + prod((M %*% alpha / alpha_0)[c(i, j)]) 
      }
    }
  }
  
  Ewij
  Ewij_emp
  
  Ewij_P
  Ewij_P_emp
  
  Ewij_shadow
  Ewij_shadow_emp
  
  sum(diag(Ewij))
  sum(diag(Ewij_shadow))
  
  

  z_bar <- sum( w_bar * z )
  tmp1 <- 0
  tmp2 <- 0
  for ( i in 1:n ) {
    tmp1 <- tmp1 + Ewij_emp[i, i] * (z[i] - z_bar)^2
    if ( i < n ) {
      for ( j in (i+1):n ) {
        tmp2 <- tmp2 + Ewij_emp[i, j] * (z[i] - z_bar) * (z[j] - z_bar)
      }
    }
  }
  ans <- tmp1 + 2 * tmp2
  ans
  var(theta)
  
  
  z_shadow_bar <- sum( w_shadow_bar * z )
  tmp1 <- 0
  tmp2 <- 0
  for ( i in 1:n ) {
    tmp1 <- tmp1 + Ewij_shadow_emp[i, i] * (z[i] - z_shadow_bar)^2
    if ( i < n ) {
      for ( j in (i+1):n ) {
        tmp2 <- tmp2 + Ewij_shadow_emp[i, j] * (z[i] - z_shadow_bar) * (z[j] - z_shadow_bar)
      }
    }
  }
  ans <- tmp1 + 2 * tmp2
  ans
  var(theta_shadow)
  
  
  
  Ewij_emp
  Ewij_shadow_emp
  
  
  

  
  
  
  # Weighted sample variance (unbiased)
  s_unb <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = alpha / n, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "unbiased" )
  S_P_unb <- cov.wt( x = matrix(z, ncol = 1), 
                     wt = alpha_prime / n, 
                     cor = FALSE, 
                     center = TRUE, 
                     method = "unbiased" )
  
  ###
  apply( t( M %*% diag(rep(1, n)) ), 2, mean )
  z_star <- apply( t( M %*% diag(rep(1, n)) ), 1, function(p) { sum(p * z) } )
  S_shadow_unb <- cov.wt( x = matrix(z_star, ncol = 1), 
                          wt = w_shadow_bar,
                          cor = FALSE, 
                          center = TRUE, 
                          method = "unbiased" )
  ###
  
  sum( (theta - mean(theta))^2 ) / length(theta)
  s_unb$cov * ( sum(beta_moments) - lpNorm( alpha / n, p = 2 )^2 )
  s_unb$cov * ( mean(lp_samples^2) - lpNorm( w_bar, p = 2 )^2 )
  s_unb$cov * mean(lp_samples_rel^2)
  
  sum( (theta_P - mean(theta_P))^2 ) / length(theta_P)
  S_P_unb$cov * ( sum(beta_moments_prime) - lpNorm( w_shadow_bar, p = 2 )^2 )
  S_P_unb$cov * ( mean(lp_P^2) - lpNorm( w_shadow_bar, p = 2 )^2 ) 
  
  sum( (theta_shadow - mean(theta_shadow))^2 ) / length(theta_shadow)
  # S_shadow_unb$cov * ( mean(lp_shadow^2) - lpNorm( w_shadow_bar, p = 2 )^2 )
  # S_shadow_unb$cov * mean( lp_shadow_rel^2 )
  # S_shadow_unb$cov * ( mean(apply( shadow[idx, ], 1, function(x) { sum(x^2) } ) ) - lpNorm( w_shadow_bar, p = 2 )^2 )
  # S_shadow_unb$cov * ( mean(lp_samples^2) - lpNorm( w_bar, p = 2 )^2 )
  var(z_star) * ( sum(beta_moments) - lpNorm( w_bar, p = 2 )^2 )
  
  
  
  sum( (theta_shadow_shifted - mean(theta_shadow_shifted))^2 ) / length(theta_shadow_shifted)
  s_unb$cov * ( mean(lp_shadow_shifted^2) - lpNorm( w_bar, p = 2 )^2 )
  
  
  
  2 / (n + 1)
  range( apply( samples[idx, ], 1, function(x) { sum(x^2) } ) )
  range( apply( shadow[idx, ], 1, function(x) { sum((solve(M) %*% x)^2) } ) )
  range( apply( shadow[idx, ], 1, function(x) { sum((x)^2) } ) )
  mean( apply( shadow[idx, ], 1, function(x) { sum((x)^2) } ) )
  
  
  
  
  
  ldens <- apply( cbind(theta, theta_shadow, theta_P, theta_shadow_shifted), 2, density )
  slolz:::plot.ldensity( ldens, col = 1:length(ldens) )
  
  
  
  
  b <- var(theta_shadow) / ( sum(w_shadow_bar * (z - z_shadow_bar)^2) * ( mean(lp_shadow^2) - lpNorm( w_shadow_bar, p = 2 )^2 ) )
  b
  1 / (1 - sum(w_shadow_bar^2))
  
  
  det( abs(M) )
  determinant( abs(M), logarithm = FALSE )
  determinant( abs(M), logarithm = TRUE )
  
  
  
  
  1 / (1-det(M)) * S_P_unb$cov * ( mean(lp_shadow^2) - lpNorm( w_shadow_bar, p = 2 )^2 )
  var(theta_shadow)
  
  
  
  det(M)
  sum(diag(Ewij_emp)) - sum(w_bar^2)
  mean( apply(samples, 1, function(w) { sum( (w - w_bar)^2 ) } ) )
  
  
  sum(diag(Ewij_shadow_emp)) - sum(w_shadow_bar^2)
  mean( apply(shadow, 1, function(w) { sum( (w - w_shadow_bar)^2 ) } ) )
  
  
  
  ##########################
  
  n_sim <- 10^4
  n <- 5
  z <- setNames( c(1:n), paste0("z", 1:n) )
  alpha <- rep(1, n)
  w_bar <- alpha / n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
  shadow <- t( apply( samples, 1, function(p) { M %*% matrix(p, ncol = 1) } ) )
  w_shadow_bar <- M %*% matrix(w_bar, ncol = 1)
  z_shadow <- apply( t( M %*% diag(rep(1, n)) ), 1, function(p) { sum(p * z) } )
  z_shadow
  t(M) %*% z
  
  w_bar_emp <- apply( samples, 2, mean )
  w_shadow_bar_emp <- apply( shadow, 2, mean )
  cbind( w_shadow_bar_emp, M %*% w_bar_emp ) # same same
  
  
  2 / (n + 1)
  mean( apply( shadow, 1, function(x) { sum(x^2) } ) )
  
  
  lp2_samples <- apply( samples, 1, function(x) { sum(x^2) } )
  lp2_shadow <- apply( shadow, 1, function(x) { sum(x^2) } )
  
  ldens <- apply( cbind(lp2_samples, lp2_shadow), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  
  
  theta <- apply( samples, 1, function(w) { sum(z * w) } )
  theta_shadow <- apply( shadow, 1, function(w) { sum(z * w) } )
  
  ldens <- apply( cbind(theta, theta_shadow), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  
  
  var(z) * (2/(n+1) - 1/n); var(theta)
  var(z_shadow) * (2/(n+1) - 1/n); var(theta_shadow)
  
  
  
  
  # Geometric random walk sampling
  
  Constraints <- constraints( selection = names(z) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint( lower = setNames(rep(0, n), names(z)),
                                               upper = setNames(rep(1, n), names(z)) )
  addConstraint(Constraints) <- orderingConstraint( Names = names(z),
                                                    ordering = 1:n,
                                                    eps = 1e-04 )
  P_tmp <- as.polytope(Constraints) 
  P <- Polytope$new( A = P_tmp@A,
                     b = P_tmp@b, 
                     sense = P_tmp@sense )
  
  interior_point <- setNames( (1:n) / sum(1:n), names(z) )
  PS <- PolytopeSampler$new( polytope = P )
  PS$setCtrl( algo = "billiard",
              interior_point = interior_point,
              n_sim = n_sim * 2,
              thin = 10^2,
              jmp = 0.1,
              adaptive_jmp = FALSE )
  # PS$polytope$isInterior( x = interior_point )
  PS$transform()
  # debugonce( PS$sample )
  # debugonce( PS$polytope$sample )
  # debugonce( har )
  PS$sample()
  PS$transformBack()
  samples_rp <- PS$samples$S
  samples_rp <- samples_rp[sample(1:nrow(samples_rp))[1:n_sim], ]
  theta_rp <- apply( samples_rp, 1, function(w) { sum(w * z) } )
  lp2_samples_rp <- apply( samples_rp, 1, function(x) { sum(x^2) } ) 

  
  colors <- 1:n
  ldens <- apply( cbind(lp2_samples, lp2_shadow, lp2_samples_rp), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE, col = colors )
  
  
  ldens <- apply( cbind(theta, theta_shadow, theta_rp), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE, col = colors )
  
  
  boxplot( cbind(theta, theta_shadow, theta_rp) )
  
  
  
  var(z) * (2/(n+1) - 1/n); var(theta)
  var(z_shadow) * (2/(n+1) - 1/n); var(theta_shadow); var(theta_rp)
  
  
  # Can we use Varsi? -- looks good
  
  q_vec <- seq( from = 1e-12, to = 0.915, length.out = 10^5 )
  z_quantile <- quantile(z, q_vec)
  FUN <- function(z0) { tail(as.numeric(varsi( mu = z, b = z0)), 1) }
  tic <- Sys.time()
  p_vec <- unlist( lapply( z_quantile, FUN ) )
  (toc_varsi <- Sys.time() - tic)
  
  plot( x = z_quantile, y = p_vec )
  
  th <- 0.889
  quantile(theta, th)
  z_quantile[ which(p_vec > th)[1] ]
  
  
  
  
  z_shadow_quantile <- quantile(z_shadow, q_vec)
  FUN <- function(z0) { tail(as.numeric(varsi( mu = z_shadow, b = z0)), 1) }
  tic <- Sys.time()
  p_vec_shadow <- unlist( lapply( z_shadow_quantile, FUN ) )
  (toc_varsi <- Sys.time() - tic)
  
  plot( x = z_shadow_quantile, y = p_vec_shadow )
  
  th <- 0.889
  quantile(theta_shadow, th)
  z_shadow_quantile[ which(p_vec_shadow > th)[1] ]
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Loop over increasing n
  # --------------------------------------------------------------------------
  
  
  n_vec <- ceiling( seq(from = 3, to = 1000, length.out = 30) )
  n_boot <- 50
  B <- 10^3 + 1
  set.seed(1:10)
  df <- 6
  x <- rt( n = max(n_vec), df = df )
  z <- x^2

  Moments <- matrix( NA, nrow = length(n_vec), ncol = 4,
                     dimnames = list(paste0("n=", n_vec), 
                                     c("M2", "M2_boot", "M2_shadow", "M2_boot_shadow")) )
  # sclmat <- Moments[ ,1]
  # colnames(sclmat) <- "scl2"
  WBPN <- Moments[ ,1:2]
  EPN <- Moments[ ,1:2]
  EPNR <- Moments[ ,1:2]
 
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    w_bar <- rep( 1/n, n )
    alpha <- w_bar * n
    lNorm <- list()
    lNorm_rel <- list()
    lNorm_rel2 <- list()
    
    # Run sampler
    bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
    lP <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = B, alpha = alpha )
      lP[[j]][is.na(lP[[j]])] <- 0
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z[1:n]) } )
    }
    
    # Squared 2-norm of average weights
    w_bar_2norm2 <- sum(w_bar^2)
    
    # Expected squared 2-norm 
    lNorm <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, p = 2) } )
    expected_2norm2 <- mean( unlist( lapply( lNorm, FUN = function(x) { mean(x^2) } ) ) )
    
    # Expected p-distance to centroid (to the power of p)
    lNorm_rel <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, x_bar = w_bar, p = 2) } )
    expected_2normrel2 <-  mean( unlist( lapply( lNorm_rel, FUN = function(x) { mean(x^2) } ) ) )
    
    # Weighted sample variance (biased)
    Vz_unb <- cov.wt( x = matrix(z[1:n], ncol = 1), 
                      wt = w_bar, 
                      cor = FALSE, 
                      center = TRUE, 
                      method = "unbiased" )
    
    # Compute variance of sample mean
    M2 <- Vz_unb$cov * ( expected_2normrel2 )
    M2_boot <- mean( apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )

    # Scaling (should converge to 1 with increasing n, but from left or right?)
    WBPN[i, 1] <- w_bar_2norm2
    EPN[i, 1] <- expected_2norm2
    EPNR[i, 1] <- expected_2normrel2
    
    
    # Shadow Dirichlet
    M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
    lP_shadow <- lapply( lP, FUN = function(P) { t( apply( P, 1, function(p) { M %*% matrix(p, ncol = 1) } ) ) } )
    bb_mat_shadow <- do.call( cbind, lapply( lP_shadow, FUN = function(P) { apply( P, 1, function(p) { sum(p * z[1:n]) } ) } ) )
    w_shadow_bar <- M %*% matrix(w_bar, ncol = 1)
    
    # Squared 2-norm of average weights
    w_shadow_bar_2norm2 <- sum(w_shadow_bar^2)
    
    # Expected squared 2-norm 
    lNorm_shadow <- lapply( lP_shadow, FUN = function(P) { apply(P, 1, lpNorm, p = 2) } )
    expected_2norm2_shadow <- mean( unlist( lapply( lNorm_shadow, FUN = function(x) { mean(x^2) } ) ) )
    
    # Expected p-distance to centroid (to the power of p)
    lNorm_shadow_rel <- lapply( lP_shadow, FUN = function(P) { apply(P, 1, lpNorm, x_bar = w_bar, p = 2) } )
    expected_2normrel2_shadow <-  mean( unlist( lapply( lNorm_shadow_rel, FUN = function(x) { mean(x^2) } ) ) )
    
    # Weighted sample variance (biased)
    Vz_unb_shadow <- cov.wt( x = matrix(z[1:n], ncol = 1), 
                             wt = w_shadow_bar, 
                             cor = FALSE, 
                             center = TRUE, 
                             method = "unbiased" )
    
    # Compute variance of sample mean
    M2_shadow <- Vz_unb$cov * ( expected_2normrel2_shadow )
    M2_boot_shadow <- mean( apply( bb_mat_shadow, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
    
    WBPN[i, 2] <- w_shadow_bar_2norm2
    EPN[i, 2] <- expected_2norm2_shadow
    EPNR[i, 2] <- expected_2normrel2_shadow
    
    Moments[i, ] <- c(M2, M2_boot, M2_shadow, M2_boot_shadow)
    
  }
  
  
  plot( Moments[ ,1], Moments[ ,2] )
  plot( Moments[ ,3], Moments[ ,4] )
  
  plot( Moments[ ,1], Moments[ ,3] )
  plot( Moments[ ,2], Moments[ ,4] )
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Regularized Pmfs
  # --------------------------------------------------------------------------
  
  I <- diag(1, n)
  ones <- rep(1, n)
  theta <- 1:n / sum(1:n)
  lambda <- 0.3
  M <- (1 - lambda) * theta %*% t(ones) + lambda * I
  M  
  
  # Uniform samples
  S <- Simplex$new( d = n )
  samples <- S$runif( n_sim = n_sim )
  z_bar <- rep(1/n, n)
  
  # Dirichlet samples
  alpha <- (1:n)^2 / sum((1:n)^2) * n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  z_bar <- alpha / n
  
  # Transform samples
  shadow <- t( apply( samples, 1, function(p) { M %*% matrix(p, ncol = 1) } ) )
  z_shadow_bar <- M %*% matrix(z_bar, ncol = 1)
  cbind( z_shadow_bar, apply( shadow, 2, mean) )
  
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", cex = 1 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", cex = 3 )
  
  # Lp-norms
  lp_samples <- apply( samples, 1, lpNorm, p = 2 )
  lp_samples_rel <- apply( samples, 1, lpNorm, p = 2, x_bar = z_bar )
  lp_shadow <- apply( shadow, 1, lpNorm, p = 2 )
  lp_shadow_rel <- apply( shadow, 1, lpNorm, p = 2, x_bar = z_shadow_bar )
  
  plot( lp_samples, lp_shadow)
  plot( lp_samples_rel, lp_shadow_rel)
  
  # Expected squared lp-norms
  mean( lp_samples^2 ); mean( lp_shadow^2 )
  mean( lp_samples_rel^2 ); mean( lp_shadow_rel^2 )
  boxplot( cbind(lp_samples, lp_shadow) )
  
  
  # Acc-rej for points with fixed lp-norms
  FUN <- function(x, expected_lp_sq, th) { abs( sum(x^2) - expected_lp_sq ) < th }
  idx_samples <- apply( samples, 1, FUN, expected_lp_sq = mean(lp_samples^2), th = 1e-03 ) 
  idx_shadow <- apply( shadow, 1, FUN, expected_lp_sq = mean(lp_shadow^2), th = 1e-03 ) 
  
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1 )
  points3d( x = samples[idx_samples, 1], y = samples[idx_samples, 2], z = samples[idx_samples, 3], col = 1, size = 3 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 2 )
  points3d( x = shadow[idx_shadow, 1], y = shadow[idx_shadow, 2], z = shadow[idx_shadow, 3], col = 1, size = 3 )
  
  
  

  # Compute statistics    
  theta <- apply( samples, 1, function(p) { sum( p * z ) } )
  theta_shadow <- apply( shadow, 1, function(p) { sum( p * z ) } )

  boxplot( cbind(theta, theta_shadow) )
  
  
  # Weighted sample variance (unbiased)
  s_unb <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = w_bar, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "unbiased" )
  S_shadow_unb <- cov.wt( x = matrix(z, ncol = 1), 
                         wt = w_shadow_bar, 
                         cor = FALSE, 
                         center = TRUE, 
                         method = "unbiased" )
  
  sum( (theta - mean(theta))^2 ) / length(theta)
  as.numeric( s_unb$cov * ( mean(lp_samples^2) - sum(w_bar^2) ) )
  as.numeric( s_unb$cov * mean(lp_samples_rel^2) )

  
  sum( (theta_shadow - mean(theta_shadow))^2 ) / length(theta_shadow)
  S_shadow_unb$cov * ( mean(lp_shadow^2) - sum(w_shadow_bar^2) )
  S_shadow_unb$cov * mean( lp_shadow_rel^2 )
  
 
  
  ldens <- apply( cbind(theta, theta_shadow), 2, density )
  slolz:::plot.ldensity( ldens, col = 1:length(ldens) )
  
  
  
  a1 <- apply( samples, 1, function(x) { sum( (x - z_bar)^2 ) } )
  a2 <- apply( samples, 1, function(x) { sum(x^2) } ) - sum(z_bar^2)
  a3 <- apply( samples, 1, function(x) { sum(x^2) } ) + sum(z_bar^2) - 2 * apply( samples, 1, function(x) { sum( x * z_bar ) } )
  cbind(a1, a2, a3)

  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ACCEPTANCE-REJECTION OF GRIDNODES WITH PROBABILITY RESCALING
  # --------------------------------------------------------------------------
  
  n <- 3
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1111)
  df <- 4
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  # Partitions
  comb( n = 2*n-1, k = n )
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # Reorder results by decreasing total probability
  ordering <- rev( order( prob_tot) ) 
  parts_mat <- parts_mat[ ,ordering]
  draws <- draws[ ,ordering]
  prob <- prob[ordering]
  prob_tot <- prob_tot[ordering]
  
  # 2-norm of partitions
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  
  # p-norm of partitions
  tmp <- lapply( 1:M, function(m) { apply( parts_mat / n, 2, lpNorm, p = m ) } )
  lp_parts_mat <- do.call( cbind, tmp )
  lp_parts_mat
  
  # Entropy
  entropy <- function( p, eps = 1e-5, exponential = TRUE, base = exp(1) )
  {
    N <- -sum(p * log(1 + (p - 1) * (p > eps), base = base))
    if (isTRUE(exponential)) {
      N <- exp(N)
    }
    return(N)
  }
  ent_parts <- apply( parts_mat / n, 2, entropy )
  ent_parts <- apply( parts_mat / n, 2, entropy, base = 2 )
  
  
  # Permutations
    lP <- list()
  lTheta <- list()
  for( j in 1:ncol(parts_mat) ) {
    
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    lP[[j]] <- perm_unique
    theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
    lTheta[[j]] <- theta
  }
  gridnodes_all <- do.call( rbind, lP )
  Grid <- gridnodes_all / n
  dim(gridnodes_all)
  comb( n = 2*n-1, k = n )
  
  # gridnodes_all
  # Grid

  # Node activation  probabilities
  lNAP <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  prob_k <- unlist(lNAP)
  sum(prob_k)
  
  # p-norm of grid nodes
  w_bar <- rep(1/n, n)
  m <- 2
  lNorm <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = m) } )
  lNorm2wbar <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = m, x_bar = w_bar) } )

  
  lp_mat <- do.call( cbind, lapply( lNorm, FUN = function(x) { unlist( lapply(x, FUN = unique) ) } ) )
  lp_mat
  
  # lp_mat_rel <- do.call( cbind, lapply( lNorm2wbar, FUN = function(x) { unlist(lapply(x, FUN = unique)) } ) )
  lp_mat_rel <- do.call( cbind, lapply( lNorm2wbar, FUN = function(x) { unlist(lapply(x, FUN = function(w) { w[1] } ) ) } ) )
  lp_mat_rel
  
  
  
  # # Sampling (Bootstrap)
  # wmat <- matrix(0, nrow = B, ncol = n)
  # for ( i in 1:B ) {
  #   idx <- sample( 1:n, replace = TRUE )
  #   tbl <- table(idx)
  #   wmat[i, as.numeric(names(tbl))] <- tbl / n
  # }
  # theta_boot <- apply( wmat, 1, function(w) { sum(w * z) } )
  # m2_theta_boot <- sum( (theta_boot - mean(theta_boot))^2 ) / length(theta_boot)
  # m2_theta_boot
  
  # Classical bootstrap by hand
  idx <- numeric(B)
  for ( i in 1:B ) {
    idx[i] <- sample( 1:nrow(Grid), size = 1, replace = TRUE, prob = prob_k )
  }
  theta_boot <- apply( Grid[idx, ], 1, function(w) { sum(w * z) } )
  m2_theta_boot <- sum( (theta_boot - mean(theta_boot))^2 ) / length(theta_boot)
  
  m2_theta_boot
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = 2/n - 1/n^2 )
  
  
  # # Sample mean
  # z_bar <- mean(z)
  # m2_hat <- var(z)
  # m2_hat_ml <- m2_hat * (n-1) / n
  # 
  # # Theoretical variance of z_bar
  # m2 / n          # using population variance
  # m2_hat / n      # using sample variace
  # m2_hat_ml / n   # using biased sample variance
  # 
  # # Bootstrap (classical)
  # n_boot <- 50
  # B <- 10^4 + 1
  # boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  # tic <- Sys.time()
  # for ( j in 1:n_boot ) {
  #   Boot <- boot::boot( data = z,
  #                       statistic = meanBoot,
  #                       R = B )
  #   boot_mat[ ,j] <- Boot$t
  # }
  # (toc_boot <- Sys.time() - tic)
  # mean( apply(boot_mat, 2, var) )
  
  
  
  mean(z)
  mean(unlist(lTheta))  # same same
  sum( nap * unlist(lTheta) ) # same same
  
 


  
  
  # Exact m2 using expectation under the distribution over the Grid
  # \mathds{E}_G (\theta - \mathds{E}_G(\theta))^2
  b <- z - mean(z)
  Ans <- matrix(0, nrow = nrow(Grid), ncol = 2)
  for ( k in 1:nrow(Grid) ) {
    g <- as.numeric(Grid[k, ])
    Ans[k, 1] <- sum( g^2 * b^2 )
    tmp <- 0
    for ( i in 1:(n-1) ) {
      for ( j in (i+1):n ) {
        tmp <- tmp + g[i] * g[j] * b[i] * b[j] * 2
      }
    }
    Ans[k, 2] <- tmp
  }
  ans <- apply( Ans, 1, sum )
  
  apply( Ans, 2, function(x) { sum( prob_k * x ) } )
  sum(  apply( Ans, 2, function(x) { sum( prob_k * x ) } ) )
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = 2/n - 1/n^2 ) # same same
  A <- b %*% t(b)
  tmp <- apply( Grid, 1, function(g) { t(g) %*% A %*% g } )
  sum( prob_k * tmp ) # same same
  m2_theta_boot  # close enough
  
  
  
  
  require(GPO)
  require(garcholz)
  Names <- aviationAlphabet(1:n)
  colnames(A) <- Names
  rownames(A) <- Names
  Constraints <- constraints(Names)
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Cons_ls) <- boxConstraint( name = "LongOnly" )
  GPS <- gps( Data = NULL,
              Covariance = A,
              Constraints = Constraints )
  GPO <- gpo(GPS)
  w_opt <- getWeights(GPO)
  Cons_ls <- Constraints
  addConstraint(Cons_ls) <- boxConstraint( name = "LongShort" )
  GPS <- gps( Data = NULL,
              Covariance = A,
              Constraints = Cons_ls )
  GPO <- gpo(GPS)
  w_opt_ls <- getWeights(GPO)
  
  cbind(w_bar, w_opt, w_opt_ls)
  
  
  t(w_bar) %*% A %*% w_bar # == zero
  t(w_opt) %*% A %*% w_opt
  t(w_opt_ls) %*% A %*% w_opt_ls
  w_tmp <- 1:n / sum(1:n)
  
  plot( apply( Grid, 1, function(g) { t(g) %*% A %*% g } ) )
  barplot( c(t(w_bar) %*% A %*% w_bar, t(w_opt) %*% A %*% w_opt) )
  
  
  S <- rdirichlet( n = 10^4, alpha = rep(1, n) )
  E <- Ellipsoid$new( centre = rep(0, n),
                      shape = make.positive.definite(A),
                      rhs = 5e-17 )
  ES <- EllipsoidSampler$new( ellipsoid = E )
  ES$setCtrl( n_sim = 10^4,
              algo = "sphere",
              b_pushy = FALSE )
  ES$sample()
  samples_E <- ES$samples$S
  
  plot3d( x = S[ ,1], y = S[ ,2], z = S[ ,3] )
  points3d( x = samples_E[ ,1], y = samples_E[ ,2], z = samples_E[ ,3], col = "grey" )
  # // Ellipsoid is degenerate ( a disc in n = 3 )
  
  
  
  # Can we use expected squared weights? --> Yes!
  wsq_bar <- apply( Grid^2, 2, function(x) { sum( prob_k * x ) } )
  sum( wsq_bar * b^2 ); sum( prob_k * Ans[ ,1] )   # same same ????????????????????

  # Expected product w_i w_j
  W <- t( apply( Grid, 1, function(g) { (g %*% t(g))[upper.tri(g %*% t(b))] } ) )
  wiwj_bar <- apply( W, 2, function(x) { sum( prob_k * x ) } )
  wisq_bar <- apply( Grid^2, 2, function(x) { sum( prob_k * x ) } )   
  wiwj_bar_mat <- matrix(NA, nrow = n, ncol = n, byrow = TRUE)
  wiwj_bar_mat[upper.tri(wiwj_bar_mat)] <- wiwj_bar
  mat <- t(wiwj_bar_mat)
  mat[upper.tri(mat)] <- wiwj_bar
  wiwj_bar_mat <- t(mat)
  wiwj_bar_mat
  tmp <- 0
  for ( i in 1:(n-1) ) {
    for ( j in (i+1):n ) {
      tmp <- tmp + wiwj_bar_mat[i, j] * b[i] * b[j] * 2
    }
  }
  tmp; mean( Ans[ ,2] )   # same same
  sum(wsq_bar * b^2) + tmp;  sum( prob_k * ans )  # almost same
 
  
  
  
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                    wt = w_bar,
                    cor = FALSE,
                    center = TRUE,
                    method = "unbiased" )
  expected_2norm2 <- sum( prob_k * apply( Grid, 1, lpNorm, p = 2 )^2 )
  Vz_unb$cov * ( expected_2norm2 - sum(w_bar^2) )
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = 2/n - 1/n^2 ) # same same
  
  expected_2norm2_rel <- sum( prob_k * apply( Grid, 1, lpNorm, x_bar = w_bar, p = 2 )^2 )
  Vz_unb$cov * expected_2norm2_rel
 

  

  
  
  
    
  
  # --------------------------------------------------------------------------
  # Constrained domain
  # --------------------------------------------------------------------------
  
  # idx_excl <- c(4, 5, 6, 10)
  # w_bar_G <- apply( Grid[-idx_excl, ], 2, mean )
  # prob_k_G <- prob_k
  # prob_k_G[idx_excl] <- 0
  # prob_k_G <- prob_k_G / sum(prob_k_G)
  
  # idx_incl <- which(Grid[ ,ncol(Grid)] >= 1 - 3/ncol(Grid))
  # idx_incl <- which( apply( Grid, 1, function(x) { max(x) == 1 } ) )
  # idx_incl <- which( apply( Grid, 1, function(x) { max(x) >= 0.5 } ) )
  idx_incl <- c(4, 6, 10)
  w_bar_G <- apply( Grid[idx_incl, ], 2, mean )
  prob_k_G <- prob_k * 0
  prob_k_G[idx_incl] <- prob_k[idx_incl]
  prob_k_G <- prob_k_G / sum(prob_k_G)
  
  cbind( prob_k, prob_k_G )
  
  sum( prob_k * unlist(lTheta) )
  sum( prob_k_G * unlist(lTheta) )
  
  
  # Classical bootstrap by hand
  B <- 10^5
  idx <- numeric(B)
  for ( i in 1:B ) {
    idx[i] <- sample( 1:nrow(Grid), size = 1, replace = TRUE, prob = prob_k_G )
  }
  theta_boot_G <- apply( Grid[idx, ], 1, function(w) { sum(w * z) } )
  m2_theta_boot_G <- sum( (theta_boot_G - mean(theta_boot_G))^2 ) / length(theta_boot_G)
  
  
  b <- z - sum(w_bar_G * z)
  Ans <- matrix(0, nrow = nrow(Grid), ncol = 2)
  for ( k in 1:nrow(Grid) ) {
    g <- as.numeric(Grid[k, ])
    Ans[k, 1] <- sum( g^2 * b^2 )
    tmp <- 0
    for ( i in 1:(n-1) ) {
      for ( j in (i+1):n ) {
        tmp <- tmp + g[i] * g[j] * b[i] * b[j] * 2
      }
    }
    Ans[k, 2] <- tmp
  }
  ans <- apply( Ans, 1, sum )
  
  
  
  tmp <- setNames( unlist( lapply( sort(unique(idx)), FUN = function(x) { sum(idx == x) } ) ),
                   sort(unique(idx)) )
  tmp
  cbind( tmp / sum(tmp), prob_k_G[prob_k_G > 0] )
  apply( cbind( tmp / sum(tmp), prob_k_G[prob_k_G > 0] ), 2, mean ) 
  
  
  apply( Ans, 2, function(x) { sum( prob_k_G * x ) } )
  sum(  apply( Ans, 2, function(x) { sum( prob_k_G * x ) } ) ); m2_theta_boot_G  # close enough
  
  A <- b %*% t(b)
  tmp <- apply( Grid, 1, function(g) { t(g) %*% A %*% g } )
  sum( prob_k_G * tmp ) # same same
  
  
  
  
  
  
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                    wt = w_bar_G,
                    cor = FALSE,
                    center = TRUE,
                    method = "unbiased" )
  l2sq_rel <- apply( Grid, 1, lpNorm, x_bar = w_bar_G, p = 2 )^2
  # l2sq_rel <- apply( Grid, 1, function(x) { sum( (x - w_bar_G)^2 ) } )
  expected_2norm2_rel_G <- sum( prob_k_G * l2sq_rel )
  Vz_unb$cov * expected_2norm2_rel_G                 # can be very different  ?????
  sum(  apply( Ans, 2, function(x) { sum( prob_k_G * x ) } ) ); m2_theta_boot_G
  
  
  
  
    
  
  
  
   