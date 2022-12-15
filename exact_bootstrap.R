  
  
  ############################################################################
  ### EXACT BOOTSTRAP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     23.12.2020
  # First version:    23.12.2020
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
  
  
  
  # Normal population, sample standard deviation
  # Varsi vs. Dirichlet
  # Compare frustum_of_simples (volesti) vs. varsi (RP)
  # Posterior of Beta random variable
  # Is the exact bootstrap distribution a stepwise approximation to the 
  # exact Bayesian bootstrap distribution?
 
  
  
  # --------------------------------------------------------------------------
  # Normal population, sample standard deviation
  # --------------------------------------------------------------------------
  
  seed <- 1234
  n <- 3
  B <- 10^3
  mu_true <- 0
  sd_true <- 1
  set.seed( seed )
  x <- rnorm(n, mu_true, sd_true)
  z <- (x - mean(x))^2
  
  # Analytic
  r <- rchisq( 10^5, df = n - 1 )
  r_scl <- r * sd_true^2 / (n - 1)
  
  
  # Ordinary bootstrap
  Boot <- boot( data = z,
                statistic = meanBoot,
                R = B )
  pdf_boot <- density(Boot$t)
  cdf_boot <- cumsum( pdf_boot$y / sum(pdf_boot$y) )
  
  
  # Exact ordinary bootstrap
  tbl <- rbind(c(3, 0, 0),
               c(0, 3, 0),
               c(0, 0, 3),
               c(2, 1, 0),
               c(1, 2, 0),
               c(0, 2, 1),
               c(0, 1, 2),
               c(2, 0, 1),
               c(1, 0, 2),
               c(1, 1, 1))
  gridpoints <- t( apply(tbl, 1, function(x) { x / sum(x) }) )
  prob <- c(1, 1, 1, 3, 3, 3, 3, 3, 3, 6) * (1/3)^3
  n_gp <- comb( n = 2*n-1, k = n )
  
  theta <- apply(gridpoints, 1, function(p) { sum(p * z) } )
  
  plot( x = theta, y = prob, type = "h", lwd = 3)
  
  # cdf
  ordering <- order(theta)
  
  plot( x = theta[ordering], y = cumsum(prob[ordering]) )
  
  
  
  # parts_mat <- partitions::parts(n)
  # draws <- gridDraws(parts_mat)
  # prob <- draws2prob(draws)
  # prob_scl <- prob / sum(prob)
  # 
  # tic <- Sys.time()
  # wmat <- matrix(0, nrow = B, ncol = n)
  # for ( i in 1:B ) {
  #   idx <- sample(1:n, replace = TRUE)
  #   tbl <- table(idx)
  #   wmat[i, as.numeric(names(tbl))] <- tbl / n
  # }
  # ( toc_boot <- Sys.time() - tic )
  # 
  # 
  # # Generate table of unique grid nodes
  # lwmat <- list()
  # for ( i in 1:nrow(wmat) ) {
  #   lwmat[[i]] <- wmat[i, ]
  # }
  # tbl_list <- lwmat[ !duplicated(lapply(lwmat, sort)) ]
  # gridpoints <- matrix( unlist(tbl_list), ncol = n, byrow = TRUE )
  # 
  # # Check if there are 2d-1 \choose d possible gridpoints
  # nrow(gridpoints)
  # comb( n = 2*n-1, k = n )
  # 
  # 
  # theta_values_gp <- apply(gridpoints, 1, function(x) { t(x) %*% z })
 
  
  
  # Varsi
  
  # q_vec <- seq( from = 0.00001, to = 0.999999, length.out = 10^3 )
  # qz <- quantile(z, q_vec)
  # qz <- seq(from = min(z), to = max(z), length.out = 10^3)
  # qz <- qz^3
  qz <- pdf_boot$x
  # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) } 
  p_vec <- unlist( lapply( qz, FUN ) )
  
  # Derive density
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(qz[idx] - qz[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  sum(d_vec_std * qz[-c(1, length(qz))])
  
  
  plot( x = qz, y = p_vec )
  lines( x = pdf_boot$x, y = cdf_boot, col = 2 )
  
  plot( x = qz[idx-1], y = d_vec )
  lines( pdf_boot, col = 2 )
  lines( density(r_scl), col = 3 )
  
  mean(r_scl); Boot$t0
  
  P <- rdirichlet( n = B, alpha = rep(1/n, n) )
  z_hat <- as.numeric(P %*% z)
  ent <- entropy( rep(1/n, n), exponential = TRUE)
  P2 <- rdirichlet( n = B, alpha = rep(1/n, n) * ent )
  z_hat2 <- as.numeric(P2 %*% z)
  
  lines(density(z_hat), col = 4)
  lines(density(z_hat2), col = 5)
  lines(density(z), col = "orange")
  
  
  # Check if the Bayesian bootstrap via Dirichlet sampling converges to 
  # the exact solution with Varsi's algorithm.
  
  B_vec <- 10^(1:6)
  lzhat <- list()
  for ( i in seq(along = B_vec) ) {
    # P <- rdirichlet( n = B_vec[i], alpha = rep(1/n, n) * n )
    U <- matrix( runif( n * B_vec[i] ), nrow = B_vec[i] )
    P <- t( apply( U, 1, function(x) { log(x) / sum(log(x)) } ) )
    lzhat[[i]] <- as.numeric(P %*% z)
  }
  
  ldens <- lapply( lzhat, FUN = density )
  names(ldens) <- paste0("^", 1:length(ldens))
  
  plot.ldensity( ldens, fillin = FALSE )
  lines( x = qz[idx-1], y = d_vec )
  lines( pdf_boot, col = 4 )
  lines( density( unlist(lzhat)), col = 5 )
  # abline( v = z, col = "grey" )
  
  
  mean(r_scl); Boot$t0; mean(z); sum(d_vec_std * qz[-c(1, length(qz))])
  unlist( lapply( lzhat, FUN = mean ) )
  
  
  
  lcdf <- lapply( ldens, FUN = function(x) { list(x = x$x, y = cumsum(x$y / sum(x$y) ) ) })
  plot.ldensity( lcdf, fillin = FALSE )
  lines( x = qz, y = p_vec )
  lines( x = pdf_boot$x, y = cdf_boot, col = 4 )
  
  
  
  plot( x = theta, y = prob  / max(prob), type = "h", lwd = 3)
  lines( x = qz[idx-1], y = d_vec / max(d_vec) )
  lines( pdf_boot, col = 4 )
  lines( ldens[[length(ldens)]], col = 5 )
  
  
  
  # Compute CDF manually
  
  samples <- sort(lzhat[[3]])
  # tmp <- unlist( lapply( qz, FUN = function(x) { length(which(samples <= x)) } ) )
  qvec <- qz
  pvec <- qvec * 0
  for ( i in seq(along = qvec) ) {
    idx <- which(samples <= qvec[i] )
    pvec[i] <- length(idx)
  }
  pvec <- pvec / max(pvec)
  
  # Derive density
  idx <- 3:length(pvec)
  dvec <- (pvec[idx] - pvec[idx-2]) / abs(qvec[idx] - qvec[idx-2])
  dvec_std <- dvec / sum(dvec)
  
  
  
  plot( qvec, pvec )
  lines( qvec, p_vec )
  
  plot( qvec[-c(1, length(qvec))], dvec_std )
  lines( qvec[-c(1, length(qvec))], d_vec_std )
  
  
  
  
 
  
  
  
  
  # Colorcoding
  
  require(rgl)
 
  P <- rdirichlet( n = 10^5, alpha = rep(1/n, n) * n )
  color <- colors <- fBasics::divPalette(n = nrow(P), "Spectral")
  val <- apply(P, 1, function(p) { p %*% z } )
  colors <- color[rank(val)]
  plot3d( x = P[ ,1], y = P[ ,2], z = P[ ,3], col = colors,
          xlab = "x", ylab = "y", zlab = "z")
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Varsi vs. Dirichlet
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(RP)
  
  n <- 3
  mu_true <- 0
  sd_true <- 1
  seed <- 1234
  set.seed( seed )
  x <- rnorm(n, mu_true, sd_true)
  z <- (x - mean(x))^2
  # z <- x
  
  # As in Cales, Chalkis, Emiris (2019)
  z <- c(0, 0.01, 0.015)
  
  
  
  # Varsi
  q_vec <- seq(from = min(z), to = max(z), length.out = 500)
  
  ## Using volesti
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0 ) }
  p_volesti <- unlist( lapply( q_vec, FUN ) )
  
  ## Using RP
  FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1 ) }
  p_rp <- unlist( lapply( q_vec, FUN ) )
  # p_vec2 <- q_vec * NA
  # for ( i in seq(along = q_vec) ) {
  #   # debugonce( varsi )
  #   A <- varsi( mu = z, b = q_vec[i] )
  #   p_vec2[i] <- tail(as.numeric(A), 1)
  # }
 
  # Derive density
  idx <- 3:length(p_volesti)
  d_volesti <- (p_volesti[idx] - p_volesti[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_volesti_std <- d_volesti / sum(d_volesti)
  d_rp <- (p_rp[idx] - p_rp[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_rp_std2 <- d_rp / sum(d_rp)
  
  # Dirichlet
  P <- rdirichlet( n = 10^5, alpha = rep(1/n, n) * n )
  z_hat <- as.numeric( P %*% z )
  samples <- sort(z_hat)
  p_ds <- q_vec * 0
  for ( i in seq(along = q_vec) ) {
    idx <- which(samples <= q_vec[i] )
    p_ds[i] <- length(idx)
  }
  p_ds <- p_ds / max(p_ds)
  
  # Derive density
  idx <- 3:length(p_ds)
  d_ds <- (p_ds[idx] - p_ds[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_ds_std <- d_ds / sum(d_ds)
  
  
  # Plot CDF
  plot( q_vec, p_ds, type = "l" )
  lines( q_vec, p_volesti, col = 2 )
  lines( q_vec, p_rp, col = 3 )
  lines( x = theta[ordering], y = cumsum(prob[ordering]), type = "o", col = 4 )
  abline( h = 0.5, col = "grey" )
  abline( v = q_vec[ which(p_ds >= 0.5)[1] ], col = "grey" )
  abline( v = q_vec[ which(p_volesti >= 0.5)[1] ], col = 2 )
  
  
  # Plot PDF
  plot( q_vec[-c(1, length(q_vec))], dvec_std, type = "l" )
  lines( q_vec[-c(1, length(q_vec))], d_vec_std, col = 2 )
  lines( q_vec[-c(1, length(q_vec))], d_vec_std2, col = 3 )
  
  
  
  # --------------------------------------------------------------------------
  # Compare frustum_of_simples (volesti) vs. varsi (RP)
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(RP)
  
  n <- 3
  set.seed(1234)
  x <- rnorm(n)
  z <- (x - mean(x))^2
  q_vec <- seq(from = min(z), to = max(z), length.out = 500)
  
  # volesti
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0 ) }
  p_vec <- unlist( lapply( q_vec, FUN ) )
  
  frustum_of_simplex( a = z, z0 = qvec[2] )
  
  # RP
  p_vec2 <- q_vec * NA
  for ( i in seq(along = q_vec) ) {
    # debugonce( varsi )
    A <- varsi( mu = z, b = q_vec[i] )
    p_vec2[i] <- tail(as.numeric(A), 1)
  }
  
  # Derivatives
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  d_vec2 <- (p_vec2[idx] - p_vec2[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec_std2 <- d_vec2 / sum(d_vec2)
  
  
  # Plots
  
  head( cbind(p_vec, p_vec2) )
  
  plot( x = q_vec, y = p_vec, type = "l" )
  lines( x = q_vec, y = p_vec2, col = 2 )
  
  # Plot PDF
  plot( q_vec[-c(1, length(q_vec))], d_vec_std, type = "l" )
  lines( q_vec[-c(1, length(q_vec))], d_vec_std2, col = 2 )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Posterior of Beta random variable
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(RP)
  require(slolz)
  
  n <- 10^3
  set.seed(1234)
  z <- rbeta( n = n, shape1 = 3, shape2 = 1 )
  n1 <- sum(z)
  mu <- rbeta( n = 10^6, shape1 = n1, shape2 = n - n1 )
  mu2 <- rbeta( n = 10^6, shape1 = 1 + n1, shape2 = 1 + n - n1 )
  
  
  # x <- rbeta( n = 10^6, shape1 = 1, shape2 = 1 ) # = uniform
  # plot(density(x))
  
  
  plot( densFUN(z) )
  lines(densFUN(mu))
  abline( v = mean(z) )

  
  
  # Varsi
  # q_vec <- seq(from = min(z), to = max(z), length.out = 500)
  q_vec <- seq( from = min(mu), to = max(mu), length.out = 300 )
  FUN <- function(z0) { tail( as.numeric(varsi( mu = z, b = z0 )), 1) }
  # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0 ) }
  p_vec <- unlist( lapply( q_vec, FUN ) )  
  
  plot( x = q_vec, y = p_vec )
  abline( v = mean(z) )
  
  # Derivative
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  
  # Dirichlet sampling
  P <- rdirichlet( n = 999, alpha = rep(1/n, n) * n )
  samples <- P %*% z
  
  
 
  dens <- densFUN(mu)
  plot(dens)
  lines(densFUN(mu2), col = 3)
  lines(densFUN(samples, from = q_vec[2], to = q_vec[length(q_vec)-1], n = length(q_vec)-2), 
        col = 4)
  lines( x = q_vec[-c(1, length(q_vec))], y = d_vec_std, col = 2 )
  abline( v = mean(z) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Is the exact bootstrap distribution a stepwise approximation to the 
  # exact Bayesian bootstrap distribution?
  # --------------------------------------------------------------------------
  
  require(volesti)
  require(RP)
  require(slolz)
  source("H:/Papers Cyril/PhD_Papers/Bootstrap/R/custom_functions.R")
  
  n <- 6
  B <- 10^6
  set.seed(1234)
  z <- rbeta( n = n, shape1 = 3, shape2 = 1 )
 
  
  # Number of gridpoints
  parts_mat <- partitions::parts(n)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_scl <- prob / sum(prob)

  tic <- Sys.time()
  wmat <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample(1:n, replace = TRUE)
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )


  # Generate table of unique grid nodes
  lwmat <- list()
  for ( i in 1:nrow(wmat) ) {
    lwmat[[i]] <- wmat[i, ]
  }
  tbl_list <- lwmat[ !duplicated(lapply(lwmat, sort)) ]
  gridpoints <- matrix( unlist(tbl_list), ncol = n, byrow = TRUE )

  # Check if there are 2d-1 \choose d possible gridpoints
  nrow(gridpoints)
  comb( n = 2*n-1, k = n )


  theta_values_gp <- apply(gridpoints, 1, function(x) { t(x) %*% z })
  
  
  
  
  

  
  

  
  
  
  
  
  