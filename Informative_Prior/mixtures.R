  
  
  ############################################################################
  ### MIXTURES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     06.08.2022
  # First version:    06.08.2022
  # --------------------------------------------------------------------------
  
  
  
  
  
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
  # Data
  # --------------------------------------------------------------------------
  
  n <- 3
  n_boot <- 10
  n_sim <- 10^3 + 1
  
  set.seed(1111)
  df <- 4
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  # --------------------------------------------------------------------------
  # MIXTURE OF DIRICHLET
  # --------------------------------------------------------------------------
  
  set.seed(1111)
  n <- 30
  n_sim <- 10^5
  alpha_vec_1 <- 1:n / sum(1:n) * n
  alpha_vec_2 <- rev(alpha_vec_1)
  # alpha_vec_2 <- runif(n)
  # alpha_vec_2 <- alpha_vec_2 / sum(alpha_vec_2) * n
  alpha_vec <- (alpha_vec_1 + alpha_vec_2) / 2
  alpha_vec  
  
  samples_1 <- rdirichlet( n = n_sim, alpha = alpha_vec_1 )
  samples_2 <- rdirichlet( n = n_sim, alpha = alpha_vec_2 )
  samples <- rdirichlet( n = n_sim, alpha = alpha_vec )
  samples_conv <- (samples_1 + samples_2) / 2
  
  mean( apply( samples_1, 1, function(x) { sum(x^2) } ) )
  mean( apply( samples_2, 1, function(x) { sum(x^2) } ) )
  mean( apply( samples, 1, function(x) { sum(x^2) } ) )
  mean( apply( samples_conv, 1, function(x) { sum(x^2) } ) )
  mean( apply( rbind(samples_1, samples_2), 1, function(x) { sum(x^2) } ) )
  2 / (n + 1)
  
  
  
  theta_1 <- apply( samples_1, 1, function(w) { sum(w * z) } )
  theta_2 <- apply( samples_2, 1, function(w) { sum(w * z) } )
  theta <- apply( samples, 1, function(w) { sum(w * z) } )
  theta_conv <- apply( samples_conv, 1, function(w) { sum(w * z) } )
  
  mean(theta_1)
  mean(theta_2)
  mean(theta)
  
  var(theta_1)
  var(theta_2)
  var(theta)
  var(theta_conv)
  var(c(theta_1, theta_2))    
      
      
      
  exp2norm2_1 <- 0
  exp2norm2_2 <- 0
  exp2norm2 <- 0
  for ( i in 1:n ) {
    exp2norm2_1 <- exp2norm2_1 + betaMoments( a = alpha_vec_1[i], b = sum(alpha_vec_1) - alpha_vec_1[i], m = 2 )
    exp2norm2_2 <- exp2norm2_2 + betaMoments( a = alpha_vec_2[i], b = sum(alpha_vec_2) - alpha_vec_2[i], m = 2 )
    exp2norm2 <- exp2norm2 + betaMoments( a = alpha_vec[i], b = sum(alpha_vec) - alpha_vec[i], m = 2 )
  }
  m2_1 <- M2Exact( z = z, w_bar = alpha_vec_1 / n, exp2norm2 = exp2norm2_1 )
  m2_2 <- M2Exact( z = z, w_bar = alpha_vec_2 / n, exp2norm2 = exp2norm2_2 )
  m2 <- M2Exact( z = z, w_bar = alpha_vec_2 / n, exp2norm2 = exp2norm2 )

  m2_1
  m2_2
  
  mean(c(m2_1, m2_2))
  mean(c(var(theta_1), var(theta_2)))
  var(theta)
  m2
  var(theta_conv)
  
  
  
  
  
  
  
  
  ########################################
  
  
  # --------------------------------------------------------------------------
  # Run sampler (classical bootstrap), n out of n
  # --------------------------------------------------------------------------
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_grid <- boot_mat
  lP_boot <- list()
  lNorm_boot <- list()
  prob <- rep(1/n, n)
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample( 1:n, replace = TRUE, prob = prob )
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    }
    lP_boot[[j]] <- gridnodes
    lNorm_boot[[j]] <- apply( lP_boot[[j]], 1, lpNorm, p = 2 )
    boot_mat_grid[ ,j] <- gridnodes %*% z
  }
  
  
  mean( apply( boot_mat, 2, var ) )
  mean( apply( boot_mat_grid, 2, var ) )
  mean( unlist( lapply( lNorm_boot, function(x) { mean(x^2)} ) ) ); 2/n - 1/n^2
  
  
  