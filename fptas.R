  
  
  ############################################################################
  ### EXACT BAYESIAN BOOTSTRAP AS FPTAS TO EXACT ORDINARY BOOTSTRAP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     27.02.2021
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(numDeriv)
  require(partitions)
  require(boot)
  require(volesti)
  require(RP)
  require(slolz)
  source("H:/Papers Cyril/PhD_Papers/Bootstrap/R/custom_functions.R")
  
  
  
  

  # NORMAL POPULATION, SAMPLE STANDARD DEVIATION


  # --------------------------------------------------------------------------
  # NORMAL POPULATION, SAMPLE STANDARD DEVIATION
  # --------------------------------------------------------------------------
  
  n <- 10
  B1 <- 10^3
  B2 <- 10^7
  mu_true <- 0
  sd_true <- 1
  seed <- 1234
  set.seed( seed )
  x <- rnorm(n, mu_true, sd_true)
  z <- (x - mean(x))^2
  
  
  # Analytic
  r <- rchisq( 10^5, df = n - 1 )
  r_scl <- r * sd_true^2 / (n - 1)
  
  
  # Ordinary bootstrap
  Boot <- boot( data = z,
                statistic = meanBoot,
                R = B1 )
  pdf_boot <- density(Boot$t)
  pdf_boot$y <- pdf_boot$y / sum(pdf_boot$y)
  cdf_boot <- cumsum( pdf_boot$y / sum(pdf_boot$y) )
  
  
  # Bayesian bootstrap (Dirichlet sampling)
  P <- rdirichlet( n = B1, alpha = rep(1, n) )
  theta_bb <- apply( P, 1, function(p) { sum(p * z) } )
  
  pdf_bb <- density(theta_bb)
  pdf_bb$y <- pdf_bb$y / sum(pdf_bb$y)
  cdf_bb <- cumsum( pdf_bb$y / sum(pdf_bb$y) )
  
  
  # Exact ordinary bootstrap
  
  # Number of gridpoints
  parts_mat <- partitions::parts(n)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_scl <- prob / sum(prob)
  
  # Sample
  tic <- Sys.time()
  wmat <- matrix(0, nrow = B2, ncol = n)
  for ( i in 1:B2 ) {
    idx <- sample(1:n, replace = TRUE)
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )
  
  ## Generate table of unique gridpoints 
  tic <- Sys.time()
  lwmat <- list()
  for ( i in 1:nrow(wmat) ) {
    lwmat[[i]] <- wmat[i, ]
  }
  lGridpoints <- lwmat[ !duplicated( lwmat ) ]
  ( toc_gridpoints <- Sys.time() - tic )
  gp_unique <- do.call( rbind, lGridpoints )
  
  ## Check if there are 2d-1 \choose d possible gridpoints
  ncol(parts_mat)
  nrow(gp_unique)
  comb( n = 2*n-1, k = n )
  
  ## Generate vector of probabilitie associeted to each unique gridpoint (p^(k))
  lP <- lapply( lGridpoints, FUN = function(x) { apply(parts_mat, 2, function(p) { all(p - rev(sort(x * n)) == 0) } ) } )
  p_k <- unlist( lapply( lP, FUN = function(p) { prob[ p ] } ) )
 
  theta_k <- unlist( lapply( lGridpoints, FUN = function(x) { t(x) %*% z } ) )
  
  
  plot( x = theta_k, y = p_k, type = "h", lwd = 3)
  
  
  ## Compute exact distribution function
  theta_k_sort <- sort(theta_k)
  d_k <- p_k * NA
  for ( i in seq(theta_k) ) {
    idx <- which(theta_k <= theta_k_sort[i])
    d_k[i] <- sum(p_k[idx])
  }
  
  plot( x = theta_k_sort, y = d_k )
  
  ## Derive exact density
  idx <- 3:length(theta_k_sort)
  d_vec <- (d_k[idx] - d_k[idx-2]) / abs(theta_k_sort[idx] - theta_k_sort[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  pdf_exact <- d_vec_std
  cdf_exact <- d_k
  pdf_exact_kernel <- density(theta_k)
  pdf_exact_kernel$y <- pdf_exact_kernel$y / sum(pdf_exact_kernel$y)
  cdf_exact_kernel <- cumsum( pdf_exact_kernel$y / sum(pdf_exact_kernel$y) )
  
 
 
  
  
  # Varsi
  
  qz <- pdf_boot$x
  # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) } 
  p_vec <- unlist( lapply( qz, FUN ) )
  
  # Derive density
  idx <- 3:length(p_vec)
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(qz[idx] - qz[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  sum(d_vec_std * qz[-c(1, length(qz))])
  
  pdf_varsi <- d_vec_std
  cdf_varsi <- p_vec 
  
  
  # Numerical derivative (using numDeriv package)
  d_vec_num <- grad( func = FUN, x = qz )
  
  
  plot( x = qz, y = p_vec )
  plot( x = qz, y = d_vec_num )
  plot( x = qz[-c(1, length(qz))], y = d_vec_std )
  
  
  
  
  
  
  
  
  
  # Plot PDF
  plot( x = pdf_boot$x, y = pdf_boot$y / sum(pdf_boot$y), type = "o" )
  points( x = theta_k_sort[-c(1, length(theta_k_sort))], y = pdf_exact, col = 3, type = "h" )
  points( x = qz[-c(1, length(qz))], y = pdf_varsi, col = 2, type = "o" )
  points( x = pdf_exact_kernel$x, y = pdf_exact_kernel$y, col = 4, type = "o" )
  points( x = pdf_bb$x, y = pdf_bb$y, col = 5, type = "o" )
  
  
  # Plot CDF
  plot( x = pdf_boot$x, y = cdf_boot )
  points( x = theta_k_sort, y = d_k, col = 3 )
  points( x = qz, y = cdf_varsi, col = 2 )
  points( x = pdf_bb$x, y = cdf_bb, col = 4 )
  # points( x = pdf_exact_kernel$x, y = cdf_exact_kernel, col = 4 )
  
  
  
  # SAVE
  
  env <- new.env()
  
  
  
  
  

  
  
  
  
  
    
  
  
  
  