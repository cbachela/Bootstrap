  

  ############################################################################
  ### EXAMPLE IN BERTSIMAS AND STURT (2020), N = 81
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.05.2021
  # First version:    23.12.2020
  # --------------------------------------------------------------------------
  
  # Replicate example in Bertsimas and Sturt (2020)
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )

  
  
  # --------------------------------------------------------------------------
  # Specifications
  # --------------------------------------------------------------------------
  
  n <- 81
  B <- 10^3
  n_boot <- 10^3
  p_th <- 0.025
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, n - length(z)))
  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
 
  # parts_mat <- partitions::parts( n )
  # draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  # prob <- draws2prob( draws = draws )
  # prob_tot <- prob * draws[2, ]
  # prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  # 
  # # Reorder results by decreasing total probability
  # ordering <- rev( order( prob_tot) ) 
  # parts_mat <- parts_mat[ ,ordering]
  # draws <- draws[ ,ordering]
  # prob <- prob[ordering]
  # prob_tot <- prob_tot[ordering]
  # prob_dmult <- prob_dmult[ordering]
  # 
  # # Euclidean norm of partitions
  # lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  # 
  # 
  # head( cbind(t(draws), log(prob), log(prob_dmult), prob_tot) )
 
  
  
    
  env <- readRDS( file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  n <- env$n
  B <- env$B
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  prob_tot <- env$prob_tot
  lp_parts <- env$lp_parts
 
  
  
  
  
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
  
  
  
  
  # Save
  env$boot_mat <- boot_mat
  env$bb_mat <- bb_mat
  saveRDS( env, file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  
  
  
  
  
  
  
  
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
  z0
  
  
  # Derivative
  idx <- 3:length(p_vec)
  # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(z_vec[idx] - z_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  plot( x = z_vec[-c(1, length(z_vec))], y = d_vec_std )
  
  sum( d_vec_std * z_vec[-c(1, length(z_vec))] )
  mean(z)
  
  
  
  
  # # --------------------------------------------------------------------------
  # # Sample from high-probability partitions
  # # --------------------------------------------------------------------------
  # 
  # idx <- which(cumsum(prob_tot) >= 0.9)[1]
  # idx
  # sum( draws[2, 1:idx] ) / sum(draws[2, ] )
  # plot( draws[2, 1:10^3 ] )
  # 
  # N <- 10^5
  # theta_mat <- matrix( NA, nrow = N, ncol = idx ) 
  # P <- matrix( NA, nrow = N, ncol = nrow(parts_mat) )
  # for ( k in 1:idx ) {
  #   for ( j in 1:N ) {
  #     P[j, ] <- sample( x = parts_mat[ ,k], replace = FALSE )
  #   }  
  #   theta_mat[ ,k] <- (P / n) %*% z
  # }
  # 
  # mu_hat <- apply( theta_mat, 2, mean )
  # sds_hat <- apply( theta_mat, 2, sd )
  # sk_hat <- apply( theta_mat, 2, skewness )
  # kurt_hat <- apply( theta_mat, 2, kurtosis )
  # 
  # 
  # plot( x = lp_parts[1:idx]^2 , y = sds_hat^2 ) #, ylim = c(0, var(z)) )   # !!!!!!
  # points( x = lp_parts[1:idx]^2, y = lp_parts[1:idx]^2 * var(z) - min(lp_parts)^2 * var(z), col = 2 )
  # 
  # 
  # 
  # theta_mat[1:5, 270:ncol(theta_mat)]
  # plot( sds_hat )
  # 
  # sum( prob_tot[1:idx] * sds_hat^2 )
  # 
  # 
  # lp_gridnodes <- apply( gridnodes, 1, lpNorm, p = 2 )
  # range(lp_gridnodes)
  # range(lp_parts[1:idx])
  
  
  
  
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
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$n <- n
  env$B <- B
  env$z <- z
  env$parts_mat <- parts_mat
  env$draws <- draws
  env$prob <- prob
  env$prob_tot <- prob_tot
  env$lp_parts <- lp_parts
  env$SE <- SE
  
  saveRDS( env, file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # mementify
  # --------------------------------------------------------------------------
  

  require(orthopolynom)
 
  n <- 81
  m_max <- 20
  x <- seq(from = 1010, to = 1070, by = 10)
  x <- c(x, rep(1, n - length(x)))
  x <- x / max(x)
  z <- scale(x, FALSE, FALSE)
  moments = unlist( lapply(1:m_max, function(i) { mean(z^i) } ) )
  xgrid <- seq(from = min(z), to = max(z), length.out = 100)
  res <- momentify(moments = moments, xgrid = xgrid )
  
  plot(res, main = paste("M =",m_max)) 
  lines(density(z), col = 3)
  
 
  
  require(Dynamic)
  X <- getMSCIData( universe = "bm", frqncy = "d" )
  BBQ <- BBSRC$new()
  BBQ$setCtrl()
  BBQ$data$X_level <- cumulated(X, "discrete")
  BBQ$runRobust()
  # states <- (1 + BBQ$output$states) / 2
  # states <- BBQ$output$states_fuzzy
  states <- (1 + BBQ$output$states_fuzzy) / 2
  
  m_max <- 20
  moments = unlist( lapply(1:m_max, function(i) { mean(states^i) } ) )
  res <- momentify(moments = moments)

  plot(res, main = paste("M =",m_max))
  plot( density(states) )
  
  
  
  states <- (1 + BBQ$output$states) / 2
  DS <- dirichletSampling( Y_train = states,
                           X_train = X,
                           X_eval = X[1, ],
                           sclfct = NULL,
                           scl_by_entropy = FALSE,
                           sclfct_lp = 1,
                           n_lag = 0,
                           n_sim = 10^5,
                           weights_fun = "kernel" )
  z <- DS$X_eval_1
  z <- (max(z) - z) / (max(z) - min(z))
  xgrid <- seq( from = min(z), to = max(z), length.out = 200 )
  m_max <- 10
  moments = unlist( lapply(1:m_max, function(i) { mean(z^i) } ) )
  res <- momentify( moments = moments, xgrid = xgrid )
  
  plot( res, main = paste("M =",m_max) ) 
  lines( density(z), col = 2 )
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    require(boot)
    require(volesti)
    require(slolz)
    require(RP)
    wd <- "H:/R/notingit/Bootstrap/"
    source( paste0(wd, "Source/custom_functions.R") )
    
    env <-  readRDS( file = paste0(wd, "waRehouse/example_BertsimasSturt2020_n=81.rds") )
    
    
    
  }
  
  
  