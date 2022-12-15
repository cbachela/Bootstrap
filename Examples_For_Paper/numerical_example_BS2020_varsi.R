  
  
  ############################################################################
  ### NUMERICAL EXAMPLE OF BERTSIMAS AND STURT 2020 REVISITED
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.09.2022
  # First version:    16.06.2022
  # --------------------------------------------------------------------------
  
  # Example in Bertsimas and Sturt (2020) - Varsi
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
 
  # --------------------------------------------------------------------------
  # Load data
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  n <- env$n
  z <- env$z
  B <- env$B
  p_th <- env$p_th
  n_boot <- env$n_boot
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  prob_tot <- env$prob_tot
  lp_parts <- env$lp_parts
  boot_mat <- env$boot_mat
  bb_mat <- env$bb_mat
  
  
  
  
  # --------------------------------------------------------------------------
  # Varsi, i.e., exact BB
  # --------------------------------------------------------------------------
  
  se <- mean( apply(boot_mat, 2, sd) )
  
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
  
  
  idx <- which(p_vec >= 0.025)[1]
  ( z_vec[idx] + z_vec[idx-1] ) / 2
  
  
  # Find the quantile H(0.025) iteratively
  varsiFindQuantile( z = z, th = 0.025 )
 
  
  
  # Derivative
  idx <- 3:length(p_vec)
  # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(z_vec[idx] - z_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  plot( x = z_vec[-c(1, length(z_vec))], y = d_vec_std )
  
  sum( d_vec_std * z_vec[-c(1, length(z_vec))] )
  mean(z)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  # Compare quantile values
  z0; p
  quantile( Boot$t, p_th )  
  mean( apply(boot_mat, 2, quantile, p_th) )
  # True value is 37.08
  
  
  
  # Means
  
  # Ordinary bootstrap
  mu_boot <- apply( boot_mat, 2, mean )
  
  # Exact ordinary bootstrap
  mu_parts <- mean(z)
  
  # Bayesian bootstrap (BB)
  mu_bb <- apply( bb_mat, 2, mean )
  
  # Exact Bayesian bootstrap
  z_eval <- z_vec[-c(1, length(z_vec))]
  mu_exact_bb <- sum( d_vec_std * z_eval )
  
  theta_hat <- setNames( c(mean(mu_boot), mu_parts, mean(mu_bb), mu_exact_bb),
                         c("avg_boot_mean", "exact_boot_mean", "avg_bb_mean", "exact_bb_mean") )
  theta_hat
  
  # Note: exact_bb_mean = mean(z) --> no bias
  
  
  
  
  # Standard error's
  
  # Ordinary bootstrap
  stderr_boot <- apply( boot_mat, 2, sd )
  
  
  # Exact ordinary bootstrap
  stderr_parts <- sqrt( sum( prob_tot * (lp_parts^2 * var(z) - min(lp_parts^2) * var(z)) ) )
  
  sum( prob_tot * lp_parts^2 )
  sum( prob_tot * lp_parts )
  sqrt(n)
  min(lp_parts^2)
  lpNorm( x = rep(1/n, n), p = 2 )^2
  
  tmp <- apply( (parts_mat / n)^2, 1, mean )
  tmp
  sum(tmp)
  barplot(tmp)
  
  
  
  
  
  # Bayesian bootstrap (BB)
  stderr_bb <- apply( bb_mat, 2, sd )
  
  
  
  # Exact Bayesian bootstrap
  z_eval <- z_vec[-c(1, length(z_vec))]
  sum( d_vec_std * z_eval ); mean(z)
  stderr_exact_bb <- sqrt( sum( d_vec_std * (z_eval)^2 ) - sum( d_vec_std * z_eval )^2 )
  
  
  # Compare
  
  SE <- setNames( c(mean(stderr_boot), stderr_parts, mean(stderr_bb), stderr_exact_bb),
                  c("avg_boot_se", "exact_boot_se", "avg_bb_se", "exact_bb_se") )
  SE
  
  
  
  ###
  
  range(stderr_boot); mean(stderr_boot)
  stderr_parts
  range(stderr_bb); mean(stderr_bb)
  stderr_exact_bb
  
  sqrt( lpNorm(rep(1/n, n), p = 2)^2 * var(z) )
  lp_bb <- apply( lP[[1]], 1, lpNorm, p = 2 )
  range( sqrt( lp_bb^2 * var(z) - lpNorm(rep(1/n, n), p = 2)^2 * var(z) ) )
  
  ###
  
  
  
  
  plot( x = z_vec[-c(1, length(z_vec))], y = d_vec_std )
  tmp <- rnorm( 10^4, mean(z), sd = stderr_exact_bb )
  lines( density( tmp), col = 2 )
  
  
  
  
  
  # Normal approximation and Cornish-Fisher expansion
  
  mu <- mean(z)
  sigma <- sqrt( sum( (z - mean(z))^2 ) / n^2 )
  sk <- sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2)
  mu2 <- sum( (z - mean(z))^2 ) / n
  mu4 <- sum( (z - mean(z))^4 ) / n
  ks <- (3 * mu2^2 / n^2 + (mu4 - 3 * mu2^2) / n^3) / ( sum( (z - mean(z))^2 / n / n ) )^2
  
  
  VaR.gauss( Data = matrix(1), mu = mu, sigma = sigma, p = 0.975 )
  VaR.mod( Data = matrix(1), mu = mu, sigma = sigma, skew = sk, kurt = ks, p = 0.975 )
  
  
  
  
  
  
  
  
  
  
  
  
  
  