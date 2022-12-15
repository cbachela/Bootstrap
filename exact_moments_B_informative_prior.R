  
  
  ############################################################################
  ### EXACT MOMENTS - CLASSICAL BOOTSTRAP - INFORMATIVE PRIOR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.10.2021
  # First version:    20.10.2021
  # --------------------------------------------------------------------------
  
  
  
  
  # Findings:
  
  # The first four (scaled and centred) moments of the sample mean
  # are computable in closed form, both when using uniform 
  # and when using non-uniform multinomial sampling probabilities.
  
  # In case 2 node activation probabilities are scrambled randomly. 
  # SE based on 2-norm is still correct. Not so for closed-form solutions.
  
  # Can we correct the scaling factor of the closed-form solutions?
  
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
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
  
  n <- 7
  n_boot <- 50
  B <- 10^4 + 1
  M <- 6
  
  set.seed(1111)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
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
  
  
  
  
  # --------------------------------------------------------------------------
  # Permutations
  # --------------------------------------------------------------------------
  
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
  dim(gridnodes_all)
  
  
  # Weights for each theta
  lW <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  w <- unlist(lW)
  sum(w)
  
  # p-norm of grid nodes
  lNorm2 <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = 2) } )
  lNorm3 <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = 3) } )
  lNorm4 <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = 4) } )
  lNorm5 <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = 5) } )
  lNorm6 <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = 6) } )

  # Node activation  probabilities
  param_multinomial <- rep(1/n, n)
  lNAP <- lapply( lP, FUN = function(P) { apply(P, 1, dmultinom, prob = param_multinomial, log = FALSE) } )
  lNAP
  nap <- unlist( lapply( lNAP, FUN = function(x) { x[1] } ) )
  cbind( unlist( lapply(lNAP, unique) ), prob )  # same same
  
 
  
  
  
  # --------------------------------------------------------------------------
  # Case 0: Classical bootstrap, i.e., uniform sampling with replacement
  # --------------------------------------------------------------------------
  
  param_multinomial <- rep(1/n, n)
  
  # Run sampler (classical bootstrap)
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_grid <- boot_mat
  lP_boot <- list()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample(1:n, replace = TRUE)
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    }
    lP_boot[[j]] <- gridnodes
    boot_mat_grid[ ,j] <- gridnodes %*% z
  }
 
  
  lNorm_boot <- list()
  expected_pnormp <- numeric(M)
  exact_expected_pnormp <- numeric(M)
  
  for ( m in 1:M ) {
    
    lNorm_boot[[m]] <- lapply( lP_boot, FUN = function(P) { apply( P, 1, lpNorm, p = m ) } )
    
    # Expected p-norm^p
    expected_pnormp[m] <- mean( unlist( lapply( lNorm_boot[[m]], FUN = function(x) { mean(x^m) } ) ) )
    
    # Exact expected p-norm^p (from partitions)
    exact_expected_pnormp[m] <- sum( prob_tot * lp_parts_mat[ ,m]^m )
  }
  
  
  expected_pnormp
  exact_expected_pnormp  # close enough
  
  sum( unlist(lNAP) * unlist(lNorm4)^4 )
  exact_expected_pnormp[4]
  expected_pnormp[4]
  
  
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = param_multinomial,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # p-norm of centroid
  lp_wbar <- unlist( lapply( 1:M, FUN = function(m) { lpNorm( param_multinomial, p = m ) } ) )
  
  # Compute variance of sample mean
  M2 <- Vz_ml$cov * (1 + (exact_expected_pnormp[2] - 1) / (1 - lp_wbar[2]^2) )
  Vz_ml$cov * (1 + (expected_pnormp[2] - 1) / (1 - lp_wbar[2]^2) )
  M2
  Vz_ml$cov / n
  mean( apply( boot_mat, 2, var ) )
  mean( apply( boot_mat_grid, 2, var ) )
  
  
  
  # Compute skewness of sample mean
  M3BFUN( x = z, M2 = M2, wghts = NULL )
  M3BFUN( x = z, M2 = NULL, wghts = NULL )
  mean( apply( boot_mat, 2, skewFUN ) )
  
  
  
  # Compute kurtosis of sample mean
  M4BFUN( x = z, M2 = M2, wghts = NULL )
  M4BFUN( x = z, M2 = NULL, wghts = NULL )
  mean( apply( boot_mat, 2, kurtFUN ) )
  

 
  
  theta <- apply( gridnodes_all / n, 1, function(p) { sum(p * z) } )
  
  # Exact M2
  sum( unlist(lNAP) * (theta - mean(theta))^2 )  
  M2  # same same
  Vz_ml$cov * (1 + (exact_expected_pnormp[2] - 1) / (1 - lp_wbar[2]^2) )
  
  
  var(z) * (lp_parts_mat[ ,2] - lp_wbar[2]^2)
  idx <- 15
  P <- lP[[idx]]
  lTheta[[idx]]
  mean(lTheta[[idx]])
  sum( (lTheta[[idx]] - mean(lTheta[[idx]]) )^2 ) / length(lTheta[[idx]])
  
  cbind( unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(unlist(lTheta)))^2 ) / length(theta) } ) ),
         var(z) * (lp_parts_mat[ ,2]^2 - lp_wbar[2]^2) )
  
  a <- unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(unlist(lTheta)))^2 ) / length(theta) } ) )  
  b <- (lp_parts_mat[ ,2]^2 - lp_wbar[2]^2) * var(z)
  plot( a, b, type = "o" )
  barplot( t(cbind(a, b)), beside = TRUE )
  
  
  # Exact M3
  sum( unlist(lNAP) * (theta - mean(theta))^3 ) / as.numeric(M2^(3/2))
  M3BFUN( x = z, M2 = M2, wghts = NULL )
  M3BFUN( x = z, M2 = NULL, wghts = NULL )
  mean( apply( boot_mat, 2, skewFUN ) )
  
  
  
  # Exact M4
  sum( unlist(lNAP) * (theta - mean(theta))^4 ) / as.numeric(M2^(4/2))
  M4BFUN( x = z, M2 = M2, wghts = NULL )
  M4BFUN( x = z, M2 = NULL, wghts = NULL )
  mean( apply( boot_mat, 2, kurtFUN ) )
  
  
  
  a <- unlist( lapply( lTheta, FUN = function(theta) { sum( (theta - mean(unlist(lTheta)))^4 ) / (length(theta) * M2^2) } ) )  
  b <- (lp_parts_mat[ ,4]^4 - lp_wbar[4]^4) * kurtFUN(z)
  plot( a, b ) #, type = "o" )
  
  cbind(a, b)
  
  
  
  A <- cbind(x = a, y = lp_parts_mat[ ,4]^4 )
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  plot( x = A[ ,"x"], y = A[ ,"y"]  )
  abline( a = intercept, b = slope )
  
  
  
  
  kurtFUN(x = z) * (1 + (exact_expected_pnormp[4] - 1) / (1 - lp_wbar[4]^4) )
  kurtFUN(x = z) * ( (1/n) + 1/kurtFUN(x = z) * (3 - 3 / n) )
  
  
  
  
  exact_expected_pnormp[4] * sum( (z - mean(z))^4 ) / (n * M2^2)
  exact_expected_pnormp[4] - lp_wbar[4]^4
  
  
  
  
  kurt.wt( x = unlist(lTheta),
           wt = w )
  lRHS4 <- lapply( lTheta, FUN = function(theta) { sum( (theta - mean(theta))^4 ) / length(theta) } )
  tmp <- unlist( lapply( lRHS4, FUN = function(x) { x / M2^(4/2) } ) )
  sum( tmp * prob_tot )
  
  M4 <- sum( tmp * prob_tot )
  M4
  
  
  plot( tmp, a )
  cbind( tmp, a )
  
  
  plot( prob_tot, lp_parts_mat[ ,4]^4 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 1: Asymmetric probabilities
  # --------------------------------------------------------------------------
  
  
  param_multinomial <- 1:n / sum(1:n)
  
  # Node activation probabilities
  lNAP <- lapply( lP, FUN = function(P) { apply(P, 1, dmultinom, prob = param_multinomial, log = FALSE) } )

  # Run sampler (weighted classical bootstrap)
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_grid <- boot_mat
  lP_boot <- list()
  lNorm_boot <- list()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        weights = param_multinomial,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample( 1:n, replace = TRUE, prob = param_multinomial )
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    }
    lP_boot[[j]] <- gridnodes
    lNorm_boot[[j]] <- apply( lP_boot[[j]], 1, lpNorm, p = 2 )
    boot_mat_grid[ ,j] <- gridnodes %*% z
  }
  
  # Expected squared 2-norm 
  expected_2normSq <- mean( unlist( lapply( lNorm_boot, FUN = function(x) { mean(x^2) } ) ) )
  
  # Exact expected squared 2-norm
  exact_expected_2normSq <- sum( unlist(lNAP) * unlist(lNorm)^2 )
  
  expected_2normSq
  exact_expected_2normSq  # close enough
  
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = param_multinomial,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of centre of mass
  lp_wbar <- lpNorm( param_multinomial, p = 2 )
  
  # Compute variance of sample mean
  M2 <- Vz_ml$cov * (1 + (exact_expected_2normSq - 1) / (1 - lp_wbar^2) )
  Vz_ml$cov * (1 + (expected_2normSq - 1) / (1 - lp_wbar^2) )
  M2
  Vz_ml$cov / n
  mean( apply( boot_mat, 2, var ) )
  mean( apply( boot_mat_grid, 2, var ) )
  
  
  
  # Compute skewness of sample mean
  M3BFUN( x = z, M2 = M2, wghts = param_multinomial )
  M3BFUN( x = z, M2 = NULL, wghts = param_multinomial ) # same same
  mean( apply( boot_mat, 2, skewFUN ) ) # close enough
  
  # Compute kurtosis of sample mean
  M4BFUN( x = z, M2 = M2, wghts = param_multinomial )
  M4BFUN( x = z, M2 = NULL, wghts = param_multinomial ) # same same
  mean( apply( boot_mat, 2, kurtFUN ) ) # close enough
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Case 2: Scrambled probabilities
  # --------------------------------------------------------------------------
  
  
  param_multinomial <- 1:n / sum(1:n)
  
  # Node activation probabilities
  lNAP <- lapply( lP, FUN = function(P) { apply(P, 1, dmultinom, prob = param_multinomial, log = FALSE) } )
  nap <- unlist( lNAP )
  
  # Scrambled node activation probabilities
  lNAPS <- lapply( lNAP, sample )
  nap_s <- unlist( lNAPS )
  
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_s_mat <- boot_mat
  lNorm_boot <- list()
  lNorm_boot_s <- list()
  for ( j in 1:n_boot ) {
    idx <- sample( 1:length(nap), size = B, replace = TRUE, prob = nap )
    lNorm_boot[[j]] <- apply( gridnodes_all[idx, ] / n, 1, lpNorm, p = 2 )
    boot_mat[ ,j] <- apply( gridnodes_all[idx, ] / n, 1, function(p) { p %*% z } )
    idx <- sample( 1:length(nap), size = B, replace = TRUE, prob = nap_s )
    lNorm_boot_s[[j]] <- apply( gridnodes_all[idx, ] / n, 1, lpNorm, p = 2 )
    boot_s_mat[ ,j] <- apply( gridnodes_all[idx, ] / n, 1, function(p) { p %*% z } )
  }
  
  # Center of mass
  wghts <- t(gridnodes_all / n) %*% nap
  wghts
  param_multinomial
  
  wghts_s <- t(gridnodes_all / n) %*% nap_s
  
  
  # Expected squared 2-norm 
  expected_2normSq <- mean( unlist( lapply( lNorm_boot, FUN = function(x) { mean(x^2) } ) ) )
  expected_2normSq_s <- mean( unlist( lapply( lNorm_boot_s, FUN = function(x) { mean(x^2) } ) ) )
  
  # Exact expected squared 2-norm
  exact_expected_2normSq <- sum( unlist(lNAP) * unlist(lNorm)^2 )
  exact_expected_2normSq_s <- sum( unlist(lNAPS) * unlist(lNorm)^2 )
  
  expected_2normSq
  exact_expected_2normSq  # close enough
  expected_2normSq_s
  exact_expected_2normSq_s  # close enough
  
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = param_multinomial,
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  Vz_ml_s <- cov.wt( x = matrix(z, ncol = 1), 
                     wt = wghts_s,
                     cor = FALSE, 
                     center = TRUE, 
                     method = "ML" )
  
  # 2-norm of centre of mass
  lp_wbar <- lpNorm( param_multinomial, p = 2 )
  lp_wbar_s <- lpNorm( wghts_s, p = 2 )
  
  
  # Compute variance of sample mean
  M2 <- Vz_ml$cov * (1 + (exact_expected_2normSq - 1) / (1 - lp_wbar^2) )
  Vz_ml$cov * (1 + (expected_2normSq - 1) / (1 - lp_wbar^2) )
  M2
  Vz_ml$cov / n
  mean( apply( boot_mat, 2, var ) )

  
  
  M2_s <- Vz_ml_s$cov * (1 + (exact_expected_2normSq_s - 1) / (1 - lp_wbar_s^2) )
  Vz_ml_s$cov * (1 + (expected_2normSq_s - 1) / (1 - lp_wbar_s^2) )
  M2_s
  Vz_ml_s$cov / n # different
  mean( apply( boot_s_mat, 2, var ) )
  
  
  # Compute skewness of sample mean
  M3BFUN( x = z, M2 = M2, wghts = param_multinomial )
  M3BFUN( x = z, M2 = NULL, wghts = param_multinomial ) # same same
  mean( apply( boot_mat, 2, skewFUN ) ) # close enough
  
  M3BFUN( x = z, M2 = M2_s, wghts = wghts_s )
  M3BFUN( x = z, M2 = NULL, wghts = wghts_s )
  mean( apply( boot_s_mat, 2, skewFUN ) ) 
  
  
  # Compute kurtosis of sample mean
  M4BFUN( x = z, M2 = M2, wghts = param_multinomial )
  M4BFUN( x = z, M2 = NULL, wghts = param_multinomial ) # same same
  mean( apply( boot_mat, 2, kurtFUN ) ) # close enough
  
  M4BFUN( x = z, M2 = M2, wghts = wghts_s )
  M4BFUN( x = z, M2 = NULL, wghts = wghts_s )
  mean( apply( boot_s_mat, 2, kurtFUN ) ) 
  
  
  theta <- apply( gridnodes_all / n, 1, function(p) { sum(p * z) } )
  sum( nap * (theta - mean(theta))^3 ) / as.numeric(M2^(3/2))
  
  
  tmp <- (theta - mean(theta))^3 / length(theta) / as.numeric(M2_s^(3/2))
  sum( nap_s * tmp )
  
  
  
  
  
  
  
  
  