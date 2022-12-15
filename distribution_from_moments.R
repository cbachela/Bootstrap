  
  
  ############################################################################
  ### APPROXIMATE DISTRIBUTION BY MOMENTS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.10.2021
  # First version:    22.10.2021
  # --------------------------------------------------------------------------
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  

  # install.packages("orthopolynom")
  # install.packages("PDQutils")
  require(orthopolynom)
  require(PDQutils)
  
  
  
  
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
  # 
  
  
  # Bertsimas Sturt 2020 example
  
  n <- 81
  B <- 10^6
  n_boot <- 10^2
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
  # Varsi, i.e., exact BB
  # --------------------------------------------------------------------------
  
  se <- mean( apply(boot_mat, 2, sd) )
  se
  mean( apply(bb_mat, 2, sd) )
  
  
  # q_vec <- seq( from = 0.1, to = 0.999, length.out = 10^5 )
  # z_vec <- quantile(z, q_vec)
  z_vec <- seq( from = 0, to = mean(z) + se * 10, length.out = 10^5 )
  # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) }
  tic <- Sys.time()
  p_vec <- unlist( lapply( z_vec, FUN ) )
  (toc_varsi <- Sys.time() - tic)
  
  plot( x = z_vec, y = p_vec )
  plot( x = q_vec, y = p_vec )
  
  
  idx <- which(p_vec >= t_th)[1]
  ( z_vec[idx] + z_vec[idx-1] ) / 2
  
  
  # Find the quantile H(0.025) iteratively
  
  err <- 1
  tol <- 1e-12
  z0 <- quantile(z, p_th)
  z_max <- max(z)
  i <- 1
  tic <- Sys.time()
  while( err > tol ) {
    p <- FUN( z0 = z0 )
    err <- abs(p - p_th)
    if ( p > p_th ) {
      z0 <- z0 - 1/(2^i) * z_max
    } else {
      z0 <- z0 + 1/(2^i) * z_max
    }
    i <- i + 1
  }
  (toc_varsi_direct <- Sys.time() - tic)
  i
  z0_varsi <- z0 
  z0_varsi
  
  
  # Derivative
  idx <- 3:length(p_vec)
  # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(z_vec[idx] - z_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  plot( x = z_vec[-c(1, length(z_vec))], y = d_vec_std )
  
  sum( d_vec_std * z_vec[-c(1, length(z_vec))] )
  mean(z)
  
  
  
  
  # --------------------------------------------------------------------------
  # Classical bootstrap
  # Approximate the quantile using PDQutils
  # based on the first five central moments of Angelova
  # --------------------------------------------------------------------------
  
  mu <- unlist( lapply( 1:5, FUN = function(i) { sum( (z_demeaned - mean(z_demeaned))^i ) / n } ) )
  m1 <- mean(z_demeaned)
  m2 <- sum( (z_demeaned - m1)^2 ) / n^2
  m3 <- sum( (z_demeaned - m1)^3 ) / n^3
  m4 <- 3 * mu[2]^2 / n^2 + (mu[4] - 3 * mu[2]^2) / n^3
  m5 <- 10 * mu[3] * mu[2] / n^3 + (mu[5] - 10 * mu[3] * mu[2]) / n^4
  moments_angelova <- c(m1, m2, m3, m4, m5)
  cumulants <- moment2cumulant( moms = moments_angelova )
  
  z0_approx_vec <- cumulants * NA
  for ( i in seq(along = cumulants) ) {
    z0_approx_vec[i] <- qapx_cf( p = p_th, 
                                 raw.cumulants = cumulants[1:i],
                                 support = c(-Inf, Inf), 
                                 lower.tail = TRUE, 
                                 log.p = FALSE )
  }
  z0_approx <- z0_approx_vec + attr(z_demeaned, "scaled:center")
  
  
  z0_boot <- apply( boot_mat, 2, quantile, 0.025 ) 
  range(z0_boot)
  mean(z0_boot); median(z0_boot); quantile( boot_mat, 0.025 )
  37.08
  z0_approx
  length( which(z0_boot < 25) ) / length(z0_boot)
  
  
  plot( x = 1:5, y = z0_approx_vec + attr(z_demeaned, "scaled:center"),
        ylab = "quantile estimate", pch = 20 )
  abline( h = 37.08 )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # Approximate the quantile using PDQutils
  # based on the first k central moments of Cales et al 2020
  # --------------------------------------------------------------------------
  
  m_max <- 5
  m1 <- mean(z_demeaned)
  m2 <- sum( (z_demeaned - m1)^2 ) / (n * (n + 1))
  moments <- unlist( lapply( 3:m_max, FUN = function(k) { attr( getMk( z = z_demeaned, k = k ), "B" ) } ) )
  moments <- c(m1, m2, moments)
  
  
  # # Sample from Varsi distribution
  # r <- sample( x = z_vec[-c(1, length(z_vec))], size = 10^6, 
  #              replace = TRUE, prob = d_vec_std )
  # # moments <- unlist( lapply(1:10, FUN = function(i) { mean(r^i) } ) )
  # cbind( moments, 
  #        unlist( lapply(1:10, FUN = function(i) { mean((r - mean(r))^i) } ) ) )
  
  
  raw_cumulants <- moment2cumulant(moments)
  z0_approx_bb <- qapx_cf( p = p_th, 
                           raw.cumulants = raw_cumulants, 
                           support = c(-Inf, Inf), 
                           lower.tail = TRUE, 
                           log.p = FALSE )
  z0_approx_bb <- z0_approx_bb + attr(z_demeaned, "scaled:center")
  
  
  z0_bb <- apply( bb_mat, 2, quantile, 0.025 )
  
  range(z0_bb)
  mean(z0_bb); median(z0_bb)
  z0_approx_bb
  z0_varsi
  

  
  
  
  # Compare Angelova moments to Cales moments
  cbind( c(m1, m2, m3, m4, m5), moments[1:5] )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Loop over increasing number of moments to include in qapx_cf (BB)
  # --------------------------------------------------------------------------
  
  
  # m_max <- 35
  m_max <- 20
  m1 <- mean(z_demeaned)
  m2 <- sum( (z_demeaned - m1)^2 ) / (n * (n + 1))
  moments <- unlist( lapply( 3:m_max, FUN = function(k) { attr( getMk( z = z_demeaned, k = k ), "B" ) } ) )
  moments <- c(m1, m2, moments)
  raw_cumulants <- moment2cumulant(moments)

  z0_approx_bb_vec <- moments * NA
  for ( i in seq(along = moments) ) {
    
    z0_approx_bb_vec[i] <- qapx_cf( p = p_th, 
                                    raw.cumulants = raw_cumulants[1:i], 
                                    support = c(-Inf, Inf), 
                                    lower.tail = TRUE, 
                                    log.p = FALSE )
  }
  
  
  plot( x = 1:m_max, y = z0_approx_bb_vec + attr(z_demeaned, "scaled:center"),
        ylab = "quantile estimate", pch = 20 )
  abline( h = z0_varsi )
  
  #// We see that the approximation error is already very small with just the first four moments
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap with informative prior
  # Approximate the quantile using PDQutils
  # based on the first 4 central moments of Cales et al 2020
  # --------------------------------------------------------------------------
  
  
  # Bayesian bootstrap with informative prior, sum(alpha) = n
  # wghts <- 1:n / sum(1:n)
  wghts <- n:1 / sum(1:n)
  # wghts <- rep(1/n, n)
  alpha <- wghts * n
  
  bb_mat_info <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    P <- rdirichlet( n = B, alpha = alpha )
    bb_mat_info[ ,j] <- apply( P, 1, function(p) { sum(p * z) } )
  }
  
  
  
  
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
  beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
  Vz_unb <- cov.wt( x = matrix(z, ncol = 1), 
                    wt = wghts, 
                    cor = FALSE, 
                    center = TRUE, 
                    method = "unbiased" )
  m1 <- sum( wghts * z )
  m2 <- Vz_unb$cov * ( sum(beta_moments) - sum(wghts^2) )
  m3 <-  M3BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  m4 <-  M4BBFUN( x = z, M2 = NA, wghts = wghts, scl = FALSE )
  moments <- c(0, m2, m3, m4)
  raw_cumulants <- moment2cumulant(moments)
  z0_approx_bb_info <- qapx_cf( p = p_th, 
                                raw.cumulants = raw_cumulants, 
                                support = c(-Inf, Inf), 
                                lower.tail = TRUE, 
                                log.p = FALSE )
  z0_approx_bb_info <- z0_approx_bb_info + m1
  
  z0_bb_info <- apply( bb_mat_info, 2, quantile, 0.025 )
  range(z0_bb_info)
  mean(z0_bb_info)
  z0_approx_bb_info
  
  
  
  
  
  
  
  
  
  