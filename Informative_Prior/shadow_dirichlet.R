  
  
  ############################################################################
  ### SHADOW DIRICHLET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.09.2022
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  # Bela Frigyik, Maya Gupta, Yihua Chen (2010) 
  # Shadow Dirichlet for Restricted Probability Modeling
  
  
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
  source( "H:/R/RP/R/class_Simplex.R" )
  source( "H:/R/RP/R/class_Ellipsoid.R" )
  source( "H:/R/RP/R/class_Sampler.R" )
  
  
  # --------------------------------------------------------------------------
  # Specifications
  # --------------------------------------------------------------------------
  
  n <- 3
  n_sim <- 10^6
  set.seed(1111)
  x <- rt(n = n, df = 4)
  z <- (x - mean(x))^2
  names(z) <- paste0("X", 1:n)
  w_bar <- rep(1/n, n)
  # w_bar <- 1:n / sum(1:n)
  
  
  
  
  ###
  
  # Upper bounds
  
  th <- 0.4
  x <- c(0.5, 0.3, 0.2)
                                                                                                                     
  y
  sum(y)
  apply( diag(rep(1, n)), 2, function(x) { M %*% x } )
  
  
  ###
  
  
  
  # --------------------------------------------------------------------------
  # Monotonic pmf
  # --------------------------------------------------------------------------
  
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
  
  Matrix::rankMatrix(M)
  
  
  # set.seed(1)
  # M <- t( rdirichlet(n = n, alpha = rep(1, n)) )
  # M <- t( rdirichlet(n = n, alpha = rep(1, n + 1)) )
  # Matrix::rankMatrix(M)
  

    
  # ###
  # Sigma <- cov(samples)
  # M <- t( eigen(Sigma)$vectors )
  # ###
  
  
  # Sample from Dirichlet
  alpha <- w_bar * n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  # Transform samples
  shadow <- t( apply( samples, 1, function(p) { M %*% p } ) )
  w_shadow_bar <- M %*% matrix(w_bar, ncol = 1)
  cbind( w_shadow_bar, apply( shadow, 2, mean) )
  
  # Transform Dirichlet parameter
  alpha_prime <- M %*% alpha
  cbind( alpha_prime / n, M %*% alpha/n, w_shadow_bar ) # same same
  
  # Sample from Dirichlet with transformed parameter
  # Note: contrary to shadow Dirichlet, the domain is still the entire simplex.
  P <- rdirichlet( n = nrow(samples), alpha = alpha_prime )
  
  
  
  # Compute expected 2-norms
 
  lp2_samples <- apply( samples, 1, lpNorm, p = 2 )
  lp2_samples_rel <- apply( samples, 1, lpNorm, p = 2, x_bar = w_bar )
  lp2_shadow <- apply( shadow, 1, lpNorm, p = 2 )
  lp2_shadow_rel <- apply( shadow, 1, lpNorm, p = 2, x_bar = w_shadow_bar )
  lp2_P <- apply( P, 1, lpNorm, p = 2 )
  
  mean(lp2_samples^2); mean(lp2_shadow^2); 2 / (n + 1)

  
  ldens <- apply( cbind(lp2_samples^2, lp2_shadow^2), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  
  ldens <- apply( cbind(lp2_samples_rel, lp2_shadow_rel), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )

  lpNorm( x = w_bar, p = 2 )^2
  lpNorm( x = w_shadow_bar, p = 2 )^2
  
  
  
  
  
  # Compute statistic
  theta <- apply( samples, 1, function(p) { sum( p * z ) } )
  theta_shadow <- apply( shadow, 1, function(p) { sum( p * z ) } )
  theta_P <- apply( P, 1, function(p) { sum( p * z ) } )
  
  boxplot( cbind(theta, theta_shadow, theta_P) )
  
  ldens <- apply( cbind(theta, theta_shadow), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  
  
  # Transform the basis
  w_base <- diag(rep(1, n))
  w_base_shadow <- t( M %*% w_base )  # transpose bcs. M is left-stochastic (i.e., columns are sum to one)
  
  # Transform the data
  # z_shadow <- apply( w_base_shadow, 1, function(p) { sum(p * z) } )
  z_shadow <- t(M) %*% z
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Check whether moments are correct
  # --------------------------------------------------------------------------
  
  betaFUN <- function(a, m) { betaMoments( a = a, a0 = sum(alpha), m = m ) }
  
  # M2
  beta_moments_2 <- unlist( lapply( alpha, FUN = betaFUN, m = 2 ) )
  Vz <- cov.wt( x = as.matrix(z, ncol = 1),
                wt = w_bar,
                method = "unbiased" )
  
  var(theta)
  Vz$cov * (sum(beta_moments_2) - sum(w_bar^2))
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments_2) )
  
  var(theta_shadow)
  sum( (theta_shadow - mean(theta_shadow))^2 ) / length(theta_shadow)
  M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments_2), M = M )
  
  
  # M3
  beta_moments_3 <- unlist( lapply( alpha, FUN = betaFUN, m = 3 ) )
  M3Exact( z = z, 
           w_bar = w_bar, 
           exp3norm3 = sum(beta_moments_3), 
           exp2norm2 = sum(beta_moments_2),
           M = M )
  sum( (theta_shadow - mean(theta_shadow))^3 ) / length(theta_shadow)
  
  
 
  
  
  
  
  # --------------------------------------------------------------------------
  # Geometric random walk sampling (monotonic pmf)
  # --------------------------------------------------------------------------
  
  Constraints <- constraints( selection = names(z) )
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
              jmp_adaptive = FALSE )
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
  
  
  
  
  # --------------------------------------------------------------------------
  # Can we use Varsi? -- yes
  # --------------------------------------------------------------------------
  
  th <- 0.025
  q_vec <- seq( from = 1e-12, to = 0.999, length.out = 10^5 )
  
  
  # Dirichlet
  z_quantile <- quantile(z, q_vec)
  quantile(theta, th)
  varsiFindQuantile( z = z, th = th, tol = 1e-03 )
  
  FUN <- function(z0) { tail(as.numeric(varsi( mu = z, b = z0)), 1) }
  p_vec <- unlist( lapply( z_quantile, FUN ) )
  z_quantile[ which(p_vec > th)[1] ]
  
  plot( x = z_quantile, y = p_vec )
  abline( v = z_quantile[ which(p_vec > th)[1] ] )
  abline( h = p_vec[ which(p_vec > th)[1] ] )
  
  
  
  
  # Shadow Dirichlet
  z_shadow_quantile <- quantile(z_shadow, q_vec)
  quantile(theta_shadow, th)
  varsiFindQuantile( z = z_shadow, th = th )
 
  
  FUN <- function(z0) { tail(as.numeric(varsi( mu = z_shadow, b = z0)), 1) }
  p_vec_shadow <- unlist( lapply( z_shadow_quantile, FUN ) )
  z_shadow_quantile[ which(p_vec_shadow > th)[1] ]
  p_vec_shadow[ which(p_vec_shadow > th)[1] ]
  
  plot( x = z_shadow_quantile, y = p_vec_shadow )
  
  
  
  

  
  # Compare empirical cdf of sampled theta's to Varsi
  ECDF <- ecdf(theta)
  plot( ECDF(z_quantile) )
  lines( p_vec, col = 2 )

  ECDFS <- ecdf(theta_shadow)
  plot( ECDFS(z_shadow_quantile) )
  lines( p_vec_shadow, col = 2 )
  
  
  
    
  
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  
  betaFUN <- function(a) { betaMoments( a = a, a0 = sum(alpha), m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  delta <- abs( lp2_samples^2 - sum(beta_moments) )
  idx <- which(delta < 1e-03)
  
  idx_H <- which( abs( theta - quantile(theta, 0.3) ) < 1e-04 )
  idx_H
  
  
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1, 
          box = FALSE, axes = TRUE, type = "p", xlab = "", ylab = "", zlab = "" )
  # rgl.bbox(xlen = 0, ylen = 0, zlen = 0, color = c('lightblue'))
  # points3d( x = P[ ,1], y = P[ ,2], z = P[ ,3], col = "green", size = 2 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 2 )
  points3d( x = samples[idx, 1], y = samples[idx, 2], z = samples[idx, 3], col = "darkgrey", size = 3 )
  points3d( x = shadow[idx, 1], y = shadow[idx, 2], z = shadow[idx, 3], col = "red", size = 3 )
  points3d( x = samples[idx_H, 1], y = samples[idx_H, 2], z = samples[idx_H, 3], col = "darkgrey", size = 3 )
  points3d( x = shadow[idx_H, 1], y = shadow[idx_H, 2], z = shadow[idx_H, 3], col = "red", size = 3 )
  # points3d( x = w_bar[1], y = w_bar[2], z = w_bar[3], size = 5, col = "darkgrey" )
  # points3d( x = w_shadow_bar[1], y = w_shadow_bar[2], z = w_shadow_bar[3], size = 5, col = "red" )
  # points3d( x = w_base[ ,1], y = w_base[ ,2], z = w_base[ ,3], size = 7, col = "darkgrey" )
  # points3d( x = w_base_shadow[ ,1], y = w_base_shadow[ ,2], z = w_base_shadow[ ,3], size = 7, col = "red" )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Regularized pmf
  # --------------------------------------------------------------------------
  
  
  ones <- rep(1, n)
  I <- diag(ones)
  lambda <- 0.6
  # s_vec <- (1:n)^3 / sum((1:n)^3)
  s_vec <- rep(1/n, n)
  M <- (1 - lambda) * s_vec %*% t(ones) + lambda * I
  M
  apply(M, 2, sum)
  
  
  
  
  # Sample from Dirichlet
  alpha <- s_vec * n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  # Transform samples
  shadow <- t( apply( samples, 1, function(p) { M %*% matrix(p, ncol = 1) } ) )

  # L2-norms
  lp2_samples <- apply( samples, 1, function(x) { sum(x^2) } )
  lp2_shadow <- apply( shadow, 1, function(x) { sum(x^2) } )
  mean(lp2_samples)
  mean(lp2_shadow)
    
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  delta <- abs( lp2_samples - sum(beta_moments) )
  idx <- which(delta < 1e-03)
  
  # Plot
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 2 )
  points3d( x = samples[idx, 1], y = samples[idx, 2], z = samples[idx, 3], col = "darkgrey", size = 3 )
  points3d( x = shadow[idx, 1], y = shadow[idx, 2], z = shadow[idx, 3], col = "red", size = 3 )
  
  
  
  

  
  # Compute statistic
  theta <- apply( samples, 1, function(p) { sum( p * z ) } )
  theta_shadow <- apply( shadow, 1, function(p) { sum( p * z ) } )
  theta_P <- apply( P, 1, function(p) { sum( p * z ) } )
  
  boxplot( cbind(theta, theta_shadow, theta_P) )
  ldens <- apply( cbind(theta, theta_shadow), 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  
  
  # Transform the basis
  w_base <- diag(rep(1, n))
  w_base_shadow <- t( M %*% w_base )  # transpose bcs. M is left-stochastic (i.e., columns are sum to one)
  
  # Transform the data
  z_shadow <- t(M) %*% z
  
  
  
  
  
  
  var(theta)
  var(z) * (2/(n+1) - 1/n)
  var(z) * mean(lp2_samples_rel^2)
  
  
  var(theta_shadow)
  var(z_shadow) * (2/(n+1) - 1/n)
  
  
  
  
  
  
  
  # # m > n:
  # 
  # m <- n + 1
  # samples <- rdirichlet( n = n_sim, alpha = rep(1, n) )
  # samples_m <- rdirichlet( n = n_sim, alpha = rep(1, m) )
  # M <- rdirichlet( n = n, alpha = rep(1, m) )
  # shadow <- t( apply( samples_m, 1, function(p) { M %*% p } ) )
  # 
  # plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 1 )
  # points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 2 )
  # 
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Shadow Multinomial
  # --------------------------------------------------------------------------
  
  
  prob <- rep(1/n, n)
  samples <- t( rmultinom( n = n_sim, size = n, prob = prob ) ) / n
  M <- do.call( cbind, lapply( 1:n, FUN = function(k) { c( rep(0, k-1), rep( 1/(n-k+1), (n-k+1) ) ) } ) )
  shadow <- t( apply( samples, 1, function(p) { M %*% p } ) )
  
  theta <- apply( samples, 1, function(p) { sum(p * z) } )
  theta_shadow <- apply( shadow, 1, function(p) { sum(p * z) } )
  
  var(z) * (2/n - 1/n^2 - 1/n )
  var(theta)  

    
  var(t(M) %*% z) * (2/n - 1/n^2 - 1/n )
  var(theta_shadow)  

  
  
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "black", size = 5 )
  points3d( x = shadow[ ,1], y = shadow[ ,2], z = shadow[ ,3], col = "orange", size = 6 )
  
  
  
  
    
  
  
  
  
  
  
  
  