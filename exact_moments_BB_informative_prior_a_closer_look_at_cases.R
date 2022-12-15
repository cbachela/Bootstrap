  
  
  ############################################################################
  ### EXACT MOMENTS - BAYESIAN BOOTSTRAP - INFORMATIVE PRIOR - A CLOSER LOOK
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.10.2021
  # First version:    20.10.2021
  # --------------------------------------------------------------------------
  

  # Sections:
  #
  # Case 0: Symmetric alpha with sum(alpha) = n. Loop over increasing n
  # Case 1: Asymmetric alpha, but sum(alpha) = n
  #         A closer look at Case 1: is E( (w - \bar{w})^m ) the same as for case 0 ?
  # Case 2: Symmetric alpha but sum(alpha) <> n. Loop over scalefactor
  
  
  
  
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
  
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 18
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1:10)
  df <- 5
  x <- rt( n = n, df = df )
  # z <- x^2
  z <- x
  
  
  
  
  
  # --------------------------------------------------------------------------
  # A closer look at case 0: Symmetric alpha with sum(alpha) = n
  # Loop over increasing n
  # --------------------------------------------------------------------------
  
  
  n_vec <- ceiling( seq(from = 3, to = 1000, length.out = 30) )
  n_boot <- 50
  B <- 10^3 + 1
  M <- 6
  
  set.seed(1:10)
  df <- 6
  x <- rt( n = max(n_vec), df = df )
  z <- x^2
  # z <- x
  ###
  # z <- rev(sort(z))
  ###
  
  Moments <- matrix( NA, nrow = length(n_vec), ncol = 12,
                     dimnames = list(paste0("n=", n_vec), 
                                     c("M2", "M2_boot", 
                                       "M3_a", "M3_b", "M3_boot", 
                                       "M4_a", "M4_b", "M4_boot",
                                       "mom3_boot", 
                                       "mom4_boot",
                                       "mom5_boot", 
                                       "mom6_boot")) )
  sclmat <- Moments[ ,1:5]
  colnames(sclmat) <- c("scl2", "scl3", "scl4a", "scl4b", "scl4denom")
  WBPN <- Moments[ ,rep(1, M)]
  EPN <- Moments[ ,rep(1, M)]
  EPNR <- Moments[ ,rep(1, M)]
  EPNR2 <- Moments[ ,rep(1, M)]
  EPNSQ <- Moments[ ,rep(1, M)]
  EEPN <- Moments[ ,rep(1, M)]
  colnames(EPN) <- paste0("boot_p", 1:M)
  colnames(EPNSQ) <- paste0("boot_pSq", 1:M)
  colnames(EEPN) <- paste0("exact_p", 1:M)
  
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
    
    # Prepare p-norm of centroid (to the power of p)
    w_bar_pnormp <- numeric(M)
    
    # Prepare expected p-norm (squared or to the power of p)
    expected_pnormp <- numeric(M)
    expected_pnormSq <- numeric(M)
    
    # Prepare expected p-distance to centroid (to the power of p)
    expected_pnormrelp <- numeric(M)
    expected_pnormrelp2 <- expected_pnormrelp
    
    # Prepare exact squared expected p-norm (to the power of p)
    exact_expected_pnormp <- numeric(M)
    
    for ( m in 1:M ) {
      
      # P-norm of average weights
      w_bar_pnormp[m] <- lpNorm( x = w_bar, p = m )^m
      
      # Expected p-norm (squared or to the power of p)
      lNorm[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, p = m) } )
      expected_pnormp[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^m) } ) ) )
      expected_pnormSq[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^2) } ) ) )
      
      # Expected p-distance to centroid (to the power of p)
      lNorm_rel[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, x_bar = w_bar, p = m) } )
      expected_pnormrelp[m] <-  mean( unlist( lapply( lNorm_rel[[m]], FUN = function(x) { mean(x^m) } ) ) )
      lNorm_rel2[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, function(p) { sum( (p - w_bar)^m ) } ) } )
      expected_pnormrelp2[m] <-  mean( unlist( lapply( lNorm_rel2[[m]], FUN = function(x) { mean(x) } ) ) )
      
      # Exact expected p-norm (to the power of p)
      betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
      beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
      exact_expected_pnormp[m] <- sum(beta_moments)
      
    }
    
    # w_bar_pnormp
    # expected_pnormp
    # exact_expected_pnormp
    # expected_pnormSq
    # expected_pnormrelp
    
    # Weighted sample variance (biased)
    Vz_ml <- cov.wt( x = matrix(z[1:n], ncol = 1), 
                     wt = w_bar, 
                     cor = FALSE, 
                     center = TRUE, 
                     method = "ML" )
    
    # Compute variance of sample mean
    M2 <- Vz_ml$cov * (1 + (exact_expected_pnormp[2] - 1) / (1 - 1/n) )
    M2_boot <- mean( apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
    
    # Compute skewness of sample mean
    M3_a <- M3BBFUN( x = z[1:n], M2 = M2, wghts = w_bar )
    M3_b <- M3BBFUN( x = z[1:n], M2 = NULL, wghts = NULL )
    M3_boot <- mean( apply( bb_mat, 2, skewFUN ) )
    
    # Compute kurtosis of sample mean
    M4_a <- M4BBFUN( x = z[1:n], M2 = M2, wghts = w_bar )
    M4_b <- M4BBFUN( x = z[1:n], M2 = NULL, wghts = w_bar )
    M4_boot <- mean( apply( bb_mat, 2, kurtFUN ) )
    
    # higher moments
    mom3_boot <- mean( apply( bb_mat, 2, momFUN, k = 3, scl = FALSE ) )
    mom4_boot <- mean( apply( bb_mat, 2, momFUN, k = 4, scl = FALSE ) )
    mom5_boot <- mean( apply( bb_mat, 2, momFUN, k = 5, scl = FALSE ) )
    mom6_boot <- mean( apply( bb_mat, 2, momFUN, k = 6, scl = FALSE ) )
    
    
    # Scaling (should converge to 1 with increasing n, but from left or right?)
    V1 <- sum(w_bar)
    V2 <- sum(w_bar^2)
    V3 <- sum(w_bar^3)
    V4 <- sum(w_bar^4)
    
    scl2 <- V1^2 / (V1^2 - V2)
    scl3 <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
    scl4denom <- ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
    scl4a <- (V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4))
    scl4b <- (3*V1^2 * (2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4))
    
    sclmat[i, ] <- c(scl2, scl3, scl4a, scl4b, scl4denom)
    WBPN[i, ] <- w_bar_pnormp
    EPN[i, ] <- expected_pnormp
    EPNR[i, ] <- expected_pnormrelp
    EPNR2[i, ] <- expected_pnormrelp2
    EPNSQ[i, ] <- expected_pnormSq
    EEPN[i, ] <- exact_expected_pnormp
    Moments[i, ] <- c(M2, M2_boot, 
                      M3_a, M3_b, M3_boot, 
                      M4_a, M4_b, M4_boot,
                      mom3_boot,
                      mom4_boot,
                      mom5_boot, 
                      mom6_boot)
    
  }
  
  
  # Compute almost exact relative p-norms for each n by running very large simulation
  EEPNR <- EPNR * NA
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    w_bar <- rep( 1/n, n )
    alpha <- w_bar * n
    P <- rdirichlet( n = 10^6, alpha = alpha )
    norm_rel <- unlist( lapply( 1:M, function(m) { mean( apply( P, 1, function(p) { sum( (p - w_bar)^m ) } ) ) } ) )
    EEPNR[i, ] <- norm_rel
    
  }
  headleft(EPNR)
  headleft(EPNR2)
  headleft(EEPNR)
  
  cbind( EPNR2[ ,3], EEPNR[ ,3] )
  
  
  
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    w_bar <- rep( 1/n, n )
    alpha <- w_bar * n
    P <- rdirichlet( n = 10^4, alpha = alpha )
    n4_rel <- mean( apply( P, 1, function(p) { sum( (p - w_bar)^4 ) } ) )
    n4_rel2 <- mean( apply( P, 1, function(p) { sum( p^4 - w_bar^4 - 4 * p * w_bar^3 - 4 * p^3 * w_bar + 6 * p^2 * w_bar^2 ) } ) )
    n4 <- mean( apply( P, 1, function(p) { sum(p^4) } ) )
    n3 <- mean( apply( P, 1, function(p) { sum(p^3) } ) )
    n2 <- mean( apply( P, 1, function(p) { sum(p^2) } ) )
    tmp1 <- mean( apply( P, 1, function(p) { sum( p * w_bar^3 ) } ) )
    tmp2 <- mean( apply( P, 1, function(p) { sum( p^3 * w_bar ) } ) )
    tmp3 <- mean( apply( P, 1, function(p) { sum( p^2 * w_bar^2 ) } ) )
    
    n4_rel; n4_rel2
    n4 + sum(w_bar^4) - 4 * tmp1 - 4 * tmp2 + 6 * tmp3
    
    tmp1; sum(w_bar^4)  
    n3 * sum(w_bar^2); tmp2
    
    tmp3; n2 * sum(w_bar^3)
    
    n4 - 3 * sum(w_bar^4) - 4 * n3 * sum(w_bar^2) + 6 * n2 * sum(w_bar^3)
    n4_rel  
    
    
    
  }
  
  
  
  
  
  m <- 2
  m2_ml_vec <- unlist( lapply( n_vec, FUN = function(n) { sum( (z[1:n] - mean(z[1:n]))^m ) / n } ) )
  plot( x = Moments[ ,"M2"], y = m2_ml_vec * sclmat[ ,m-1] * EPNR[ ,m] )  
  cbind( Moments[ ,"M2"], m2_ml_vec * sclmat[ ,m-1] * EPNR[ ,m] )

  
  
  
  m <- 3
  m3_ml_vec <- unlist( lapply( n_vec, FUN = function(n) { sum( (z[1:n] - mean(z[1:n]))^m ) / n } ) )
  plot( x = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2), 
        y = m3_ml_vec * sclmat[ ,m-1] * EPNR[ ,m], type = "o" )
  plot( x = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2), 
        y = m3_ml_vec * sclmat[ ,m-1] * EPNR2[ ,m], type = "o" )
  plot( x = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2), 
        y = m3_ml_vec * sclmat[ ,m-1] * EEPNR[ ,m], type = "o" )
  plot( x = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2), 
        y = m3_ml_vec * sclmat[ ,m-1] * EEPN[ ,m] )
  tmp <- cbind( m3 = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2), 
                # y = m3_ml_vec * sclmat[ ,m-1] * EPNR[ ,m] )
                y = m3_ml_vec * sclmat[ ,m-1] * EPNR2[ ,m] )
  # y = m3_ml_vec * sclmat[ ,m-1] * EEPNR[ ,m] )
  plot( log(tmp[ ,1]), log(tmp[ ,2]) )
  plot( tmp[ ,1], tmp[ ,2] )
  tmp
  ordering <- order(tmp[ ,1])
  tmp[ordering, ]
  # tmp <- tmp[ordering, ]
  barplot( t(tmp[ordering, ]), beside = TRUE )
  
  # barplot( t(tmp), beside = TRUE )
  
  
  
  
  
  
  # Checks
  cbind( Moments[ ,c("M3_a", "M3_b")], Moments[ ,"mom3_boot"] / Moments[ ,"M2"]^(3/2) )
  cbind( Moments[ ,c("M3_a", "M3_b")] * Moments[ ,"M2"]^(3/2), Moments[ ,"mom3_boot"] )
  
  FUN <- function(n)
  {
    w <- rep(1/n, n)
    ans <- 1 / (1 - 3 * sum(w^2) + 2 * sum(w^3) )
    return( ans )
  }
  cbind( sclmat[ ,2],
         unlist( lapply( n_vec, FUN = FUN ) ) )
  
  
  
  
  
  
  x_bar <- rep( 1 / 3, 3 )
  P <- rdirichlet( n = B, alpha = x_bar * length(x_bar) )
  P_scl <- t( apply( P, 1, function(p) { p - x_bar } ) ) 
  
  apply( P, 1, lpNorm, p = 2 )
  apply( P_scl, 1, lpNorm, p = 2 )
  
  require(rgl)
  plot3d( x = P[ ,1], y = P[ ,2], z = P[ ,3] )
  points3d( x = P_scl[ ,1], y = P_scl[ ,2], z = P_scl[ ,3], col = 2 )
  
  
  
  
  m <- 4
  m4 <- unlist( lapply( n_vec, FUN = function(n) { sum( (z[1:n] - mean(z[1:n]))^m ) / n } ) )
  m2 <- unlist( lapply( n_vec, FUN = function(n) { sum( (z[1:n] - mean(z[1:n]))^2 ) / n } ) )
  
  
  
  plot( x = (Moments[ ,"M4_a"] * Moments[ ,"M2"]^2)[-1],
        y = ( (sclmat[ ,"scl4a"] / sclmat[ ,"scl4denom"] * m4 - 
                 sclmat[ ,"scl4b"] / sclmat[ ,"scl4denom"] * m2^2) * 
                EPNR[ ,m] )[-1] )
  
  
  plot( x = (Moments[ ,"M4_a"] * Moments[ ,"M2"]^2)[-1],
        y = (sclmat[ ,"scl4a"] / sclmat[ ,"scl4denom"] * m4 * EPNR[ ,m] - 
               sclmat[ ,"scl4b"] / sclmat[ ,"scl4denom"] * m2^2 * EPNR[ ,m-2] )[-1] )
  
  
  
  
  
  # Exact \mathds{E}(\sum_{i=1}^n (z - \bar{z})^4 )
  exact_expected_4norm4_rel <- EEPN[ ,4] - 3*WBPN[ ,4] - 4*EEPN[ ,3]*WBPN[ ,2] + 6*EEPN[ ,2]*WBPN[ ,3]
  exact_expected_2norm2_rel <- EEPN[ ,2] - WBPN[ ,2]
  
  
  plot( x = (Moments[ ,"M4_a"] * Moments[ ,"M2"]^2)[-1],
        y = (sclmat[ ,"scl4a"] / sclmat[ ,"scl4denom"] * m4 * exact_expected_4norm4_rel - 
               sclmat[ ,"scl4b"] / sclmat[ ,"scl4denom"] * m2^2 * exact_expected_2norm2_rel )[-1] )
  
  plot( x = (Moments[ ,"M4_a"] * Moments[ ,"M2"]^2)[-1],
        y = ((sclmat[ ,"scl4a"] / sclmat[ ,"scl4denom"] * m4 - 
                sclmat[ ,"scl4b"] / sclmat[ ,"scl4denom"] * m2^2) * exact_expected_4norm4_rel )[-1] )
  
  
  
  cbind( EEPNR[ ,4],
         exact_expected_4norm4_rel )  # close enough
  plot( x = EEPNR[-1 ,4], y = exact_expected_4norm4_rel[-1] )
  
  
  
  
  
  
  plot( x = n_vec, y = Moments[ ,"M2"] )
  plot( x = n_vec, y = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2) )
  plot( x = n_vec, y = Moments[ ,"M4_a"] * Moments[ ,"M2"]^2 )
  
  plot( x = n_vec, y = m2_ml_vec )
  plot( x = n_vec, y = m3_ml_vec )

  
  plot( x = m2_ml_vec, y = Moments[ ,"M2"] )
  plot( x = m3_ml_vec, y = Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2) )
  plot( x = m4, y = Moments[ ,"M4_a"] * Moments[ ,"M2"]^2 )
  
  plot( Moments[ ,"M2"] / m2_ml_vec )
  plot( Moments[ ,"M3_a"] * Moments[ ,"M2"]^(3/2) / m3_ml_vec )
  plot( Moments[ ,"M4_a"] * Moments[ ,"M2"]^2 / m4 )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # A closer look at Case 1
  # --------------------------------------------------------------------------
  
  
  # n_boot <- 50
  # B <- 10^4 + 1
  n_boot <- 5
  B <- 10^2 + 1
  M <- 6
  
  # Dirichlet parameter (symmetric)
  wghts <- rep(1/n, n)
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
  # Generate set of asymmetric Dirichlet parameters
  S <- Simplex$new( d = n )
  alpha_mat <- S$runif( 100 )
  w_bar_mat <- alpha_mat * NA
  
  # Loop over starting points
  Moments <- matrix( NA, nrow = nrow(alpha_mat), ncol = 12,
                     dimnames = list(paste0("alpha", 1:nrow(alpha_mat)), 
                                     c("M2", "M2_boot", 
                                       "M3_a", "M3_b", "M3_boot", 
                                       "M4_a", "M4_b", "M4_boot",
                                       "mom3_boot", 
                                       "mom4_boot",
                                       "mom5_boot", 
                                       "mom6_boot")) )
  WBPN <- Moments[ ,rep(1, M)]
  EPN <- Moments[ ,rep(1, M)]
  EPNR <- Moments[ ,rep(1, M)]
  EPNSQ <- Moments[ ,rep(1, M)]
  EEPN <- Moments[ ,rep(1, M)]
  colnames(EPN) <- paste0("boot_p", 1:M)
  colnames(EPNSQ) <- paste0("boot_pSq", 1:M)
  colnames(EEPN) <- paste0("exact_p", 1:M)
  
  for ( i in 1:nrow(alpha_mat) ) {
 
    alpha <- alpha_mat[i, ] * n
    w_bar <- alpha / sum(alpha)
    w_bar_mat[i, ] <- w_bar
    lNorm <- list()
    lNorm_rel <- list()
    
    # Run sampler
    bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
    lP <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = B, alpha = alpha )
      lP[[j]][is.na(lP[[j]])] <- 0
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    # Prepare p-norm of centroid (to the power of p)
    w_bar_pnormp <- numeric(M)
    
    # Prepare expected p-norm (squared or to the power of p)
    expected_pnormp <- numeric(M)
    expected_pnormSq <- numeric(M)
    
    # Prepare expected p-distance to centroid (to the power of p)
    expected_pnormrelp <- numeric(M)

    # Prepare exact squared expected p-norm (to the power of p)
    exact_expected_pnormp <- numeric(M)
    
    for ( m in 1:M ) {
      
      # P-norm of average weights
      w_bar_pnormp[m] <- lpNorm( x = w_bar, p = m )^m
      
      # Expected p-norm (squared or to the power of p)
      lNorm[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, p = m) } )
      expected_pnormp[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^m) } ) ) )
      expected_pnormSq[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^2) } ) ) )
      
      # Expected p-distance to centroid (to the power of p)
      lNorm_rel[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, x_bar = w_bar, p = m) } )
      expected_pnormrelp[m] <-  mean( unlist( lapply( lNorm_rel[[m]], FUN = function(x) { mean(x^m) } ) ) )
      
      # Exact expected p-norm (to the power of p)
      betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
      beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
      exact_expected_pnormp[m] <- sum(beta_moments)
      
    }
    
    # w_bar_pnormp
    # expected_pnormp
    # exact_expected_pnormp
    # expected_pnormSq
    # expected_pnormrelp

    # Compute variance of sample mean
    M2 <- Vz_ml$cov * (1 + (exact_expected_pnormp[2] - 1) / (1 - lp_wbar^2) )
    M2_boot <- mean( apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
    
    # Compute skewness of sample mean
    M3_a <- M3BBFUN( x = z, M2 = M2, wghts = wghts )
    M3_b <- M3BBFUN( x = z, M2 = NULL, wghts = NULL )
    M3_boot <- mean( apply( bb_mat, 2, skewFUN ) )
  
    # Compute kurtosis of sample mean
    M4_a <- M4BBFUN( x = z, M2 = M2, wghts = wghts )
    M4_b <- M4BBFUN( x = z, M2 = NULL, wghts = wghts )
    M4_boot <- mean( apply( bb_mat, 2, kurtFUN ) )
    
    # higher moments
    mom3_boot <- mean( apply( bb_mat, 2, momFUN, k = 3, scl = FALSE ) )
    mom4_boot <- mean( apply( bb_mat, 2, momFUN, k = 4, scl = FALSE ) )
    mom5_boot <- mean( apply( bb_mat, 2, momFUN, k = 5, scl = FALSE ) )
    mom6_boot <- mean( apply( bb_mat, 2, momFUN, k = 6, scl = FALSE ) )
    
    
    WBPN[i, ] <- w_bar_pnormp
    EPN[i, ] <- expected_pnormp
    EPNR[i, ] <- expected_pnormrelp
    EPNSQ[i, ] <- expected_pnormSq
    EEPN[i, ] <- exact_expected_pnormp
    Moments[i, ] <- c(M2, M2_boot, 
                      M3_a, M3_b, M3_boot, 
                      M4_a, M4_b, M4_boot,
                      mom3_boot,
                      mom4_boot,
                      mom5_boot, 
                      mom6_boot)
    
  }
  
  
  
  
  
  #####################
  
 
  
  
  # Rimoldini_2013 Weighted skewness and kurtosis unbiased by sample size
  
  
  V1 <- sum(w_bar)
  V2 <- sum(w_bar^2)
  V3 <- sum(w_bar^3)
  V4 <- sum(w_bar^4)
  
 
  # M1
  sum(z) / n
 
  
  # M2
  scl <- V1^2 / (V1^2 - V2)
  1 / scl
  n / (n - 1)  # same same
  1 / (1 - lpNorm(w_bar, p = 2)^2)  # same same
  # w_tmp <- runif(n)
  # 1 / (1 - lpNorm(w_tmp / sum(w_tmp), p = 2)^2)  # not the same as with w_bar
  
  sum( (z - mean(z))^2 ) / n * scl
  var(z)
  
  Moments[i, "M2"]
  sum( w_bar * (z - mean(z))^2 ) / (n + 1)
  sum( w_bar * (z - mean(z))^2 ) * (1 + (EEPN[i, 2] - 1) / (1 - lpNorm(w_bar, p = 2)^2) )
  sum( w_bar * (z - mean(z))^2 ) * EPNR[i, 2] / (1 - lpNorm(w_bar, p = 2)^2)
  sum( w_bar * (z - mean(z))^2 ) * scl * EPNR[i, 2]
  sum( w_bar * (z - mean(z))^2 ) * EEPN[i, 3]
  

    
  
  
  
  
  
  
  # M3
  scl3 <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
  1 / scl3
  1 / (n^2 / ((n-1) * (n-2)))  # same same
  
  sum( (z - mean(z))^3 ) / n * scl3
  sum( (z - mean(z))^3 ) / n * (n^2 / ((n-1) * (n-2)))
  
  
  Moments[1, "mom3_boot"]
  M3BBFUN( x = z, M2 = 1, scl = FALSE )
  sum( (z - mean(z))^3 ) / n * 2 / ((n + 1) * (n + 2))   # same same
  sum( (z - mean(z))^3 ) / n * EPNR[1, 3] / (1 - lpNorm(w_bar, p = 3)^3)
  sum( (z - mean(z))^3 ) / n * scl3 * EPNR[1, 3]
 
  
  
  
  
  
  # M4
  (V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4)) / ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
  
  (3*V1^2 * (2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4)) / ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
  
  
  
  
  m2 <- sum( (z - mean(z))^2 ) / n
  m4 <- sum( (z - mean(z))^4 ) / n
  (n * (n^2 - 2*n + 3)) / ((n - 1) * (n - 2) * (n - 3)) * m4 - (3*n * (2*n - 3)) / ( (n - 1) * (n - 2) * (n - 3) ) * m2^2
  
  
  
  
  # --------------------------------------------------------------------------
  # A closer look at Case 1: Difference between moments under case 1 and case 0
  # --------------------------------------------------------------------------
  
  betaFUN <- function(a, m = 2) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
  n <- 10^2
  set.seed(1:10)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  w_bar <- 1:n / sum(1:n)
  alpha <- w_bar * n  
  
  
  # M1
  m1 <- mean(z)
  m1_w <- sum( w_bar * z )
  
  
  # M2
  m2_b <- M2Exact( z = z, exp2norm2 = 2/n - 1/n^2 )
  m2_bb <- M2Exact( z = z, exp2norm2 = 2/(n + 1) )
  beta_moments <- unlist( lapply( w_bar * n, FUN = betaFUN ) ) 
  m2_bb_w <- M2Exact( z = z, w_bar = w_bar, exp2norm2 = sum(beta_moments) )
  m2_bb
  m2_bb_w
  m2_tilde <- sum( (z - m1)^2 ) / n
  m2_tilde_w <- sum( w_bar * (z - m1_w)^2 )
  m2_tilde_w / m2_tilde; m2_bb_w / m2_bb # same same
  
  
  # M3
  m3_bb <- M3BBFUN( x = z, M2 = 1, scl = FALSE )
  m3_bb_w <- M3BBFUN( x = z, M2 = 1, wghts = w_bar, scl = FALSE )
  m3_bb
  m3_bb_w
  m3_tilde <- sum( (z - m1)^3 ) / n
  m3_tilde_w <- sum( w_bar *  (z - m1_w)^3 )
  m3_tilde
  m3_tilde_w
  m3_tilde_w / m3_tilde; m3_bb_w / m3_bb # same same
  
  
  # M4
  m4_bb <- M4BBFUN( x = z, M2 = 1, scl = FALSE )
  m4_bb_w <- M4BBFUN( x = z, M2 = 1, wghts = w_bar, scl = FALSE )
  m4_bb
  m4_bb_w
  m4_tilde <- sum( (z - m1)^4 ) / n
  m4_tilde_w <- sum( w_bar *  (z - m1_w)^4 )
  m4_tilde
  m4_tilde_w
  (6 * m4_tilde_w + 3 * n * m2_tilde_w^2) / (6 * m4_tilde + 3 * n * m2_tilde^2); m4_bb_w / m4_bb  # same same
  
  

  # Question: Is E( (w - \bar{w})^m ) the same as for case 0 ?
  # Answer: No.???
  
  n_vec <- 3:100
  m_vec <- 2:6
  betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
  
  for ( m in sq(along = m_vec) ) {
    
    for ( i in seq(along = n_vec) ) {
      
      n <- n_vec[i]
      w_bar <- rep(1/n, n)
      alpha <- w_bar * n
      # sum ( unlist( lapply( alpha, FUN = betaFUN) ) ) - sum(w_bar^m)
      2 / (n + 1) - 1 / n
      # P <- rdirichlet( n = 10^4, alpha = rep(1, n) )
      # mean( apply( P, 1, function(x) { sum( (x - w_bar)^2 ) } ) )
      Alpha <- rdirichlet( n = 10^2, alpha = rep(1, n) ) * n
      
      for ( k in 1:nrow(Alpha) )  {
        
        alpha <- as.numeric(Alpha[k, ] )
        beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
        exact_expected_m_norm_m <- sum(beta_moments)  
        tmp <- exact_expected_m_norm_m - sum((alpha/n)^m)
        tmp / biasCorrection( w_bar = alpha/n, m = m )
        tmp * biasCorrection( w_bar = alpha/n, m = m )
        
      }
    }
  }
 
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # A closer look at Case 2: (symmetric alpha, but sum(alpha) <> n)
  # Loop over scalefactor
  # --------------------------------------------------------------------------
  
  
  # Data - Bertsimas Sturt 2020 example
  n <- 81
  z <- seq(from = 1010, to = 1070, by = 10)
  z <- c(z, rep(1, n - length(z)))
  z_demeaned <- scale(z, TRUE, FALSE)
  
  # Parameters
  n_boot <- 50
  B <- 10^4 + 1
  M <- 6
  
  
  # Dirichlet parameter (symmetric)
  wghts <- rep(1/n, n)
  
  # Weighted sample variance (biased)
  Vz_ml <- cov.wt( x = matrix(z, ncol = 1), 
                   wt = wghts, 
                   cor = FALSE, 
                   center = TRUE, 
                   method = "ML" )
  
  # 2-norm of average weights
  lp_wbar <- lpNorm( x = wghts, p = 2 )
  
  # Loop over scalefactors
  # scl_vec <- seq(from = 1, to = n * 2, by = 1 )
  # scl_vec <- c( seq(from = n, to = 6/n^2, length.out = 20), c(5:1)/n^2 )
  scl_vec <- seq(from = n*2, to = 0.01, length.out = 100) 
  Moments <- matrix( NA, nrow = length(scl_vec), ncol = 12,
                     dimnames = list(paste0("scl=", scl_vec), 
                                     c("M2", "M2_boot", 
                                       "M3_a", "M3_b", "M3_boot", 
                                       "M4_a", "M4_b", "M4_boot",
                                       "mom3_boot", 
                                       "mom4_boot",
                                       "mom5_boot", 
                                       "mom6_boot")) )
  WBPN <- Moments[ ,rep(1, M)]
  EPN <- Moments[ ,rep(1, M)]
  EPNR <- Moments[ ,rep(1, M)]
  EPNR2 <- Moments[ ,rep(1, M)]
  EPNSQ <- Moments[ ,rep(1, M)]
  EEPN <- Moments[ ,rep(1, M)]
  colnames(EPN) <- paste0("boot_p", 1:M)
  colnames(EPNSQ) <- paste0("boot_pSq", 1:M)
  colnames(EEPN) <- paste0("exact_p", 1:M)
  
  for ( i in seq(along = scl_vec) ) {
    
    alpha <- wghts * scl_vec[i]
    w_bar <- alpha / sum(alpha)
    lNorm <- list()
    lNorm_rel <- list()
    lNorm_rel2 <- list()
    
    # Run sampler
    bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
    lP <- list()
    for ( j in 1:n_boot ) {
      lP[[j]] <- rdirichlet( n = B, alpha = alpha )
      lP[[j]][is.na(lP[[j]])] <- 0
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    # Prepare p-norm of centroid (to the power of p)
    w_bar_pnormp <- numeric(M)
    
    # Prepare expected p-norm (squared or to the power of p)
    expected_pnormp <- numeric(M)
    expected_pnormSq <- numeric(M)
    
    # Prepare expected p-distance to centroid (to the power of p)
    expected_pnormrelp <- numeric(M)
    expected_pnormrelp2 <- numeric(M)
    
    # Prepare exact squared expected p-norm (to the power of p)
    exact_expected_pnormp <- numeric(M)
    
    for ( m in 1:M ) {
      
      # P-norm of average weights
      w_bar_pnormp[m] <- lpNorm( x = w_bar, p = m )^m
      
      # Expected p-norm (squared or to the power of p)
      lNorm[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, p = m) } )
      expected_pnormp[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^m) } ) ) )
      expected_pnormSq[m] <- mean( unlist( lapply( lNorm[[m]], FUN = function(x) { mean(x^2) } ) ) )
      
      # Expected p-distance to centroid (to the power of p)
      lNorm_rel[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, lpNorm, x_bar = w_bar, p = m) } )
      expected_pnormrelp[m] <- mean( unlist( lapply( lNorm_rel[[m]], FUN = function(x) { mean(x^m) } ) ) )
      lNorm_rel2[[m]] <- lapply( lP, FUN = function(P) { apply(P, 1, function(p) { sum( (p - w_bar)^m ) } ) } )
      expected_pnormrelp2[m] <-  mean( unlist( lapply( lNorm_rel2[[m]], FUN = function(x) { mean(x) } ) ) )
      
      # Exact expected p-norm (to the power of p)
      betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = m ) }
      beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
      exact_expected_pnormp[m] <- sum(beta_moments)
      
    }
    
    # w_bar_pnormp
    # expected_pnormp
    # exact_expected_pnormp
    # expected_pnormSq
    # expected_pnormrelp
    
    
    # Compute variance of sample mean
    M2 <- Vz_ml$cov * (1 + (exact_expected_pnormp[2] - 1) / (1 - lp_wbar^2) )
    M2_boot <- mean( apply( bb_mat, 2, function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
    
    # Compute skewness of sample mean
    M3_a <- M3BBFUN( x = z, M2 = M2, wghts = wghts )
    M3_b <- M3BBFUN( x = z, M2 = NULL, wghts = NULL )
    M3_boot <- mean( apply( bb_mat, 2, skewFUN ) )
    
    # Compute kurtosis of sample mean
    M4_a <- M4BBFUN( x = z, M2 = M2, wghts = wghts )
    M4_b <- M4BBFUN( x = z, M2 = NULL, wghts = wghts )
    M4_boot <- mean( apply( bb_mat, 2, kurtFUN ) )
    
    # higher moments
    mom3_boot <- mean( apply( bb_mat, 2, momFUN, k = 3, scl = FALSE ) )
    mom4_boot <- mean( apply( bb_mat, 2, momFUN, k = 4, scl = FALSE ) )
    mom5_boot <- mean( apply( bb_mat, 2, momFUN, k = 5, scl = FALSE ) )
    mom6_boot <- mean( apply( bb_mat, 2, momFUN, k = 6, scl = FALSE ) )
    
    
    WBPN[i, ] <- w_bar_pnormp
    EPN[i, ] <- expected_pnormp
    EPNR[i, ] <- expected_pnormrelp
    EPNR2[i, ] <- expected_pnormrelp2
    EPNSQ[i, ] <- expected_pnormSq
    EEPN[i, ] <- exact_expected_pnormp
    Moments[i, ] <- c(M2, M2_boot, 
                      M3_a, M3_b, M3_boot, 
                      M4_a, M4_b, M4_boot,
                      mom3_boot,
                      mom4_boot,
                      mom5_boot, 
                      mom6_boot)
    
  }
  
  
  
  
  # Rimoldini_2013 Weighted skewness and kurtosis unbiased by sample size
  
  V1 <- sum(w_bar)
  V2 <- sum(w_bar^2)
  V3 <- sum(w_bar^3)
  V4 <- sum(w_bar^4)
  
  
  # M1
  sum(z) / n
  
  
  # M2
  scl2 <- V1^2 / (V1^2 - V2)
  scl2
  n / (n - 1)  # same same
  1 / (1 - lpNorm(w_bar, p = 2)^2)  # same same
  # w_tmp <- runif(n)
  # 1 / (1 - lpNorm(w_tmp / sum(w_tmp), p = 2)^2)  # not the same as with w_bar
  
  sum( (z - mean(z))^2 ) / n * scl2
  var(z)
  
  
  Moments[1, "M2"]
  sum( (z - mean(z))^2 ) / (n * (n + 1))
  sum( (z - mean(z))^2 ) / n * (1 + (EEPN[1, 2] - 1) / (1 - lpNorm(w_bar, p = 2)^2) )  # same same
  sum( (z - mean(z))^2 ) / n * 1 / (1 - lpNorm(w_bar, p = 2)^2) * (EEPN[1, 2] - lpNorm(w_bar, p = 2)^2) # same same
  sum( (z - mean(z))^2 ) / n * EPNR[1, 2] / (1 - lpNorm(w_bar, p = 2)^2) # slightly different bcs EPNR is estimated and not exact
  sum( (z - mean(z))^2 ) / n * scl2 * EPNR[1, 2]
  
  # Using exact relative norm
  EEPNR2 <- EEPN[ ,2] - lpNorm(w_bar, p = 2)^2
  sum( (z - mean(z))^2 ) / n * scl2 * EEPNR2[1]
  sum( (z - mean(z))^2 ) / (n * (n + 1))  
  
  
  m2_ml <- sum( (z - mean(z))^2 ) / n
  m2 <- m2_ml * scl2
  # M2 <- Moments[ ,"M2"] / (sum( (z - mean(z))^2 ) / n * scl2)
  M2 <- m2 * EPNR[ ,2]
  cbind( M2, EEPNR2 )  # correct
  
  
  plot( x = scl_vec, y = M2 )
  
  
  
  
  
  # M3
  scl3 <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
  scl3
  n^2 / ((n-1) * (n-2))  # same same
  1 / (1 - 3 * sum(w_bar^2) + 2 * sum(w_bar^3) )  # same same
  
  
  sum( (z - mean(z))^3 ) / n * scl3
  sum( (z - mean(z))^3 ) / n * (n^2 / ((n-1) * (n-2)))
  
  
  exact_expected_3norm3_rel <- EEPN[ ,3] - 3*EEPN[ ,2]*WBPN[ ,2] + 2*WBPN[ ,3]
  
  Moments[1, "mom3_boot"]
  M3BBFUN( x = z, M2 = 1, scl = FALSE )
  sum( (z - mean(z))^3 ) / n * 2 / ((n + 1) * (n + 2))   # same same
  sum( (z - mean(z))^3 ) / n * scl3 * as.numeric(exact_expected_3norm3_rel[1]) # same same
  sum( (z - mean(z))^3 ) / n * scl3 * EPNR2[1, 3] # close enough
  
  plot( x = Moments[ ,"mom3_boot"],
        y = EPNR2[ ,3] )
  plot( x = Moments[ ,"mom3_boot"],
        y = sum( (z - mean(z))^3 ) / n * scl3 * EPNR2[ ,3] )
  
  
  
  tmp <- cbind( x = Moments[ ,"mom3_boot"],
                y = sum( (z - mean(z))^3 ) / n * scl3 * exact_expected_3norm3_rel ) 
  tmp  
  
  barplot( t(tmp), beside = TRUE, col = 1:2 )
  plot( tmp[ ,1] / tmp[ ,2] )
  
  
  m3_ml <- sum( (z - mean(z))^3 ) / n
  m3 <- m3_ml * scl3
  M3 <-  m3 * exact_expected_3norm3_rel
  
  plot( x = scl_vec, y = map201(val = M2, max_val = m2), type = "o" )
  points( x = scl_vec, y = map201(val = M3, max_val = m3), col = 2, type = "o" )

  #// drops very rapidly
  
  
  
  
  
  # M4
  scl4denom <- ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
  scl4a <- (V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4))
  scl4b <- (3*V1^2 * (2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4))
  
  m2 <- sum( (z - mean(z))^2 ) / n
  m4 <- sum( (z - mean(z))^4 ) / n
  (n * (n^2 - 2*n + 3)) / ((n - 1) * (n - 2) * (n - 3)) * m4 - (3*n * (2*n - 3)) / ( (n - 1) * (n - 2) * (n - 3) ) * m2^2
  m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom # same same
  
  
  Moments[1, "mom4_boot"]
  Moments[1, "M4_a"] * Moments[1 ,"M2"]^2
  M4BBFUN( x = z, M2 = 1, scl = FALSE )
  (m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom) * EPNR[1, 4]
  
  plot( x = Moments[ ,"mom4_boot"],
        y = EPNR[ ,4] )  
  plot( x = Moments[ ,"mom4_boot"],
        y = (m4 * scl4a - m2^2 * scl4b) * EPNR[ ,4] )
  
  cbind( x = Moments[ ,"mom4_boot"],
         y = (m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom) * EPNR[ ,4] )
  
  
  cbind( x = Moments[ ,"mom4_boot"],
         y = EPNR[ ,4] )  
  
  
  
  # Exact \mathds{E}(\sum_{i=1}^n (z - \bar{z})^4 )
  exact_expected_4norm4_rel <- EEPN[ ,4] - 3*WBPN[ ,4] - 4*EEPN[ ,3]*WBPN[ ,2] + 6*EEPN[ ,2]*WBPN[ ,3]
  exact_expected_2norm2_rel <- EEPN[ ,2] - WBPN[ ,2]
  
  # Compare estimated with exact \mathds{E}(\sum_{i=1}^n (z - \bar{z})^4 )
  tmp <- cbind( a = EPNR[ ,4],
                b = exact_expected_4norm4_rel )
  plot( tmp[ ,1], tmp[ ,2] )
  tmp # close enough
  
  
  
  plot( x = Moments[ ,"mom4_boot"],
        y = (m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom) * exact_expected_4norm4_rel )
  plot( x = Moments[ ,"mom4_boot"],
        y = exact_expected_4norm4_rel )
  
  
  tmp <- cbind( x = Moments[ ,"mom4_boot"],
                y = (m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom) * exact_expected_4norm4_rel,
                z = (m4 * scl4a * exact_expected_4norm4_rel / scl4denom - m2^2 * scl4b * exact_expected_2norm2_rel^2 / scl4denom)) 
  cbind( tmp, tmp[ ,1] / tmp[ ,2], tmp[ ,1] / tmp[ ,3] )
  
  
  
  barplot( t(tmp), beside = TRUE, col = 1:ncol(tmp) )
  plot( x = scl_vec, y = tmp[ ,1] / tmp[ ,2] )
  plot( x = scl_vec, y = tmp[ ,1] - tmp[ ,2] )  # mom4_boot is always larger
  plot( x = scl_vec, y = tmp[ ,1] - tmp[ ,3] )  # mom4_boot is always larger
  plot( x = scl_vec, y = tmp[ ,1] / tmp[ ,3] )
  
  
  
  M4 <- Moments[ ,"mom4_boot"]
  g <- scl4a * exact_expected_4norm4_rel * m4 / scl4denom
  cc <- (g - M4) * scl4denom / (scl4b * m2^2)
  
  cc
  plot( cc, exact_expected_4norm4_rel )
  plot( cc, exact_expected_2norm2_rel^2 )
  
  barplot(cc)
  plot( cc )
  tmp <- cbind( exact_expected_4norm4_rel, exact_expected_2norm2_rel^2 )
  tmp
  
  
  f <- M4 / (m4 * scl4a / scl4denom - m2^2 * scl4b / scl4denom)
  plot(f)
  tmp <- cbind( f, exact_expected_4norm4_rel )
  cbind( tmp, tmp[ ,1] / tmp[ ,2] )
  plot( tmp[ ,1] / tmp[ ,2] )
  plot( tmp[ ,1], tmp[ ,2] )
  plot( 1/EEPN[ ,4], tmp[ ,1] / tmp[ ,2] )
  
  
  
  
  m4_ml <- sum( (z - mean(z))^3 ) / n
  m3 <- m3_ml * scl3
  M3 <-  m3 * exact_expected_3norm3_rel
  
  plot( x = scl_vec, y = map201(val = M2, max_val = m2) )
  points( x = scl_vec, y = map201(val = M3, max_val = m3), col = 2 )
  points( x = scl_vec, y = map201(val = M4, max_val = m4), col = 3 )
  
  
  
  
  
  
  
  
  
  
  #####################
  
  plot( x = scl_vec, y = Moments[ ,"M2"], ylim = c(0, var(z)), type = "o" )
  abline( h = var(z) * (n-1) / n )
  points( x = scl_vec, y = sum( (z - mean(z))^2) / n * (1 + (EEPN[ ,2] - 1) / (1 - lpNorm(w_bar, p = 2)^2) ), type = "o", col = 2 )
  
  
  plot( x = scl_vec, y = Moments[ ,"mom3_boot"], type = "o" )
  abline( h = momFUN(z, k = 3, scl = FALSE) )
  points( x = scl_vec, y = sum( (z - mean(z))^3) / n * scl3 * exact_expected_3norm3_rel, type = "o", col = 3 )
  
  
  plot( x = scl_vec, y = Moments[ ,"mom4_boot"], type = "o" )
  abline( h = momFUN(z, k = 4, scl = FALSE) )
  points( x = scl_vec, y =  (m4 * scl4a - m2^2 * scl4b) * exact_expected_4norm4_rel, type = "o", col = 3 )
  
  
  
  
  plot( x = EEPN[ ,2], y = Moments[ ,"M2"], type = "o" )
  plot( x = EEPN[ ,2], y = Moments[ ,"M2_boot"], type = "o" )
  plot( x = EEPN[ ,3], y = Moments[ ,"mom3_boot"], type = "o" )
  plot( x = EPNR[ ,3], y = Moments[ ,"mom3_boot"], type = "o" )
  plot( x = EEPN[ ,4], y = Moments[ ,"mom4_boot"], type = "o" )
  plot( x = EPNR[ ,4], y = Moments[ ,"mom4_boot"], type = "o" )
  plot( x = EEPN[ ,4] * sum( (z - mean(z))^4) / n, y = Moments[ ,"mom4_boot"], type = "o" )
  plot( x = EEPN[ ,4]^(1/4), y = Moments[ ,"mom4_boot"], type = "o" )
  plot( x = EEPN[ ,5], y = Moments[ ,"mom5_boot"], type = "o" )
  plot( x = EEPN[ ,6], y = Moments[ ,"mom6_boot"], type = "o" )
  
  
  a <- Moments[ ,"M2_boot"]
  b <- sum( (z - mean(z))^2) / n * (1 + (EEPN[ ,2] - 1) / (1 - lpNorm(w_bar, p = 2)^2) )
  
  a <- Moments[ ,"mom3_boot"]
  b <- sum( (z - mean(z))^3) / n * (1 + (EEPN[ ,3] - 1) / (1 - lpNorm(w_bar, p = 3)^3) )
  b <- sum( (z - mean(z))^3) / n * EEPN[ ,3]
  
  a <- Moments[ ,"mom4_boot"]
  b <-  sum( (z - mean(z))^4) / n * (1 + (EEPN[ ,4] - 1) / (1 - lpNorm(w_bar, p = 4)^4) )
  b <- sum( (z - mean(z))^4) / n * EEPN[ ,4]
  
  cbind(a, b)
  plot( a / b )
  cbind(a, b)[1:5, ]
  M3BBFUN(z, M2 = 1, scl = FALSE)
  
  A <- Moments[ ,"mom4_boot"] / EEPN[ ,4]^(1/4)
  plot(A)
  head(A)
  plot( EEPN[ ,4]^(1/4), A )
  
  
  A <- Moments[ ,"mom4_boot"] - EEPN[ ,4] # * sum( (z - mean(z))^4) / n
  A <- Moments[ ,"mom4_boot"] - EEPN[ ,4] * sum( (z - mean(z))^4) / n
  plot( A )
  head(A)
  plot( EEPN[ ,4], A )
  plot( EEPN[1:9, 4], A[1:9] )
  plot( Moments[ ,"mom4_boot"], EEPN[ ,4] * sum( (z - mean(z))^4) / n )
  
  
  ###
  
  Moments2 <- Moments
  Moments2[1, "mom3_boot"] <- M3BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  Moments2[1, "mom4_boot"] <- M4BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  
  fx <- Moments2[ ,"mom4_boot"] / ( sum( (z - mean(z))^4 ) / n )
  plot( x = EEPN[ ,4], y = fx, type = "o" )
  points( x = EPNR[ ,4], y = fx, col = 2, type = "o")
  points( x = fx, y = fx, col = 2, type = "o" )
  plot(fx)
  delta <- fx - EEPN[ ,4]
  plot( delta )
  plot( EEPN[ ,4], delta )
  
  tmp <- cbind( a = Moments2[ ,"mom4_boot"], 
                b = sum( (z - mean(z))^4 ) / n * EEPN[ ,4],
                delta = delta )
  head(tmp)
  
  EEPN[1, 4]
  fx[1]
  EPNR[1, 4]
  
  
  
  
  
  
  w <- runif(n)
  w <- w / sum(w)
  lpNorm(x = w, x_bar = w_bar, p = 2)^2
  lpNorm(x = w, p = 2)^2 - lpNorm(x = w_bar, p = 2)^2
  
  lpNorm(x = w, x_bar = w_bar, p = 4)^4
  lpNorm(x = w, p = 4)^4 - lpNorm(x = w_bar, p = 4)^4
  
  
  
  
  
  
  ###
  
  
  tmp <- cbind( x = EEPN[ ,4], 
                y = Moments[ ,"mom4_boot"] ) #- EEPN[ ,4] )
  slope <- (tmp[1, "y"] - tmp[nrow(tmp), "y"]) / (tmp[1, "x"] - tmp[nrow(tmp), "x"])
  intercept <- tmp[nrow(tmp), "y"] - slope * tmp[nrow(tmp), "x"]
  plot( tmp[ ,"x"], tmp[ ,"y"] )
  abline( a = intercept, b = slope )
  
  ans <- cbind( tmp[ ,"x"], intercept + slope * tmp[ ,"x"] + tmp[ ,"x"] )
  plot( ans )
  cbind( Moments[ ,"mom4_boot"], ans )
  plot( Moments[1:9 ,"mom4_boot"] )
  points( ans[1:9 ,2], col = 2 )
  
  
  
  
  
  
  
  y <- sum( (z - mean(z))^4) / n
  tmp <- cbind( Moments[ ,"mom4_boot"],
                EEPN[ ,4] * y,
                EEPN[ ,2]^2 * y,
                (EEPN[ ,4] - EEPN[ ,2]^2) * y )
  plot( as.timeSeries(tmp[1:9, ]), plot.type = "single" )
  
  head(tmp)
  
  
  
  
  k <- 4
  Moments2 <- Moments
  Moments2[1, "mom3_boot"] <- M3BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  Moments2[1, "mom4_boot"] <- M4BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  y <- Moments2[ ,paste0("mom", k, "_boot")]
  # y <- Moments[ ,"M2"]
  # y <- Moments[ ,"M2_boot"]
  x <- EEPN[ ,k]
  
  plot( x = x, y = y )
  
  A <- cbind( x = x, y = y )
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  b <- 1 - w_bar_pnormp[k]
  mz <- sum( (z - mean(z))^k) / n
  b <- (A[1, "x"] - 1) / (A[1, "y"] / mz - 1)
  b
  1 - lpNorm( w_bar, p = k )^k
  slope
  intercept
  tmp <- cbind( A, 
                A[ ,"x"] * slope + intercept, 
                mz + mz * A[ ,"x"] / b - mz / b )
  tmp
  barplot( t(tmp[1:10, -1]), beside = TRUE )
  
  
  sum( (z - mean(z))^2 ) / (length(z) * (length(z) + 1))
  Moments[1, "M2"]
  
  M4BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  Moments[1, "mom4_boot"]
  (sum( (z - mean(z))^4) / n) * (1 + (EEPN[1, 4] - 1) / (1 - w_bar_pnormp[4]) )
  
  
  M3BBFUN( x = z, M2 = Moments[1, "M2"], scl = FALSE )
  Moments[1, "mom3_boot"]
  
  
  
  #####################
  
  m2 <- sum( (z - mean(z))^2 ) / ( n * (n + 1) )
  m2
  2 * var(z) / (n + 1) - var(z) / n
  (sum( (z - mean(z))^2) / n) * (1 + (EEPN[1, 2] - 1) / (1 - w_bar_pnormp[2]) )
  a <- sum( (z - mean(z))^2 ) / n
  b <- m2 / a
  b  
  
  
  m4 <- M4BBFUN(x = z, M2 = 1, scl = FALSE)
  m4
  (sum( (z - mean(z))^4) / n) * (1 + (EEPN[1, 4] - 1) / (1 - w_bar_pnormp[4]) )
  a <- sum( (z - mean(z))^4 ) / n
  b <- m4 / a
  b  
  EEPN[1, 4]
  
  
  
  
  
  #####################
  # M2
  
  plot( x = n * Moments[ ,"M2"] / sum((z - mean(z))^2), 
        y = EEPN[ ,2], ylim = c(0, max(EEPN[ ,2])), xlim = c(0, max(EEPN[ ,2])) )
  
  
  y <- Moments[ ,"M2"]
  x <- sum((z - mean(z))^2) / n * EEPN[ ,2] 
  plot( x = x, y = y )
  A <- cbind( x = x, y = y )
  A
  
  
  y <- n * Moments[ ,"M2"] / sum((z - mean(z))^2) 
  x <- EEPN[ ,2] 
  plot( x = x, y = y )
  A <- cbind( x = x, y = y )
  A
  
  plot( x = x, y = y )
  abline(a = intercept, b = slope)
  
  
  cc <- 1/n * sum((z - mean(z))^2)
  b <- cc / (1 - w_bar_pnormp[2])
  
  
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  slope
  intercept
  cbind(A, 
        A[ ,"x"] * slope + intercept, 
        as.numeric(Vz_ml$cov) * (1 + (EEPN[ ,2] - 1) / (1 - w_bar_pnormp[2]) ),
        cc - b + b * EEPN[ ,2] )
  
  
  plot( x = x, y = y )
  abline(a = intercept, b = slope)
  
  
  
  
  ################
  # M4 
  
  M4BBFUN(x = z)
  Moments[10, "M4_boot"]
  
  Moments[10, "M4_boot"] <- M4BBFUN(x = z)
  plot( x = n * Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 / sum((z - mean(z))^4), 
        y = EEPN[ ,4] ) # !
  
  x <- n * Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 / sum( (z - mean(z))^4 )
  y <- EEPN[ ,4]
  
  # x <- EEPN[ ,4] * sum( (z - mean(z))^4 ) / (n * Moments[ ,"M2"]^2)
  # y <- Moments[ ,"M4_boot"]
  
  plot( x = x, y = y )
  
  A <- cbind( x = x, 
              y = y )
  A
  
  
  cbind( A, 
         A[ ,"x"] * slope + intercept )
  
  
  
  A <- cbind( x = (EEPN[ ,4]),
              y = Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 )
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  plot( EEPN[ ,4], Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 )
  abline(a = intercept, b = slope)
  
  
  
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  slope
  intercept
  cbind( A, 
         A[ ,"x"] * slope + intercept, 
         (EEPN[ ,4] - intercept) * sum( (z - mean(z))^4 ) / (n * Moments[ ,"M2"]^2 * slope),
         Moments[ ,"M4_boot"] )
  
  plot( x = (EEPN[ ,4] - intercept) * sum( (z - mean(z))^4 ) / (n * Moments[ ,"M2"]^2 * slope),
        y = Moments[ ,"M4_boot"] )
  
  
  plot( x = x, y = y )
  abline(a = intercept, b = slope)
  
  
  
  
  
  
  plot( x = as.numeric(Vz_ml$cov) * (1 + (EEPN[ ,2] - 1) / (1 - w_bar_pnormp[2]) ), 
        y = Moments[ ,"M2"] )
  cbind( x = as.numeric(Vz_ml$cov) * (1 + (EEPN[ ,2] - 1) / (1 - w_bar_pnormp[2]) ), 
         y = Moments[ ,"M2"] )
  
  
  
  
  plot( x = (1 + (EEPN[ ,3] - 1) / (1 - w_bar_pnormp[3]) ), 
        y = Moments[ ,"M3_boot"] )
  
  plot( x = 10^3 * Moments[ ,"M3_boot"] * Moments[ ,"M2"]^(3/2) / sum((z - mean(z))^3) )
  plot( x = 10^3 * Moments[ ,"M3_boot"] * Moments[ ,"M2"]^(3/2) / sum((z - mean(z))^3), y = exp(EEPN[ ,2]), type = "o" )
  plot( x = 10^3 * Moments[ ,"M3_boot"] * Moments[ ,"M2"]^(3/2) / sum((z - mean(z))^3), 
        y = EEPN[ ,3], type = "o" )
  
  plot( x = (10^3 * Moments[ ,"M3_boot"]),
        y = EEPN[ ,3] * sum((z - mean(z))^3) / Moments[ ,"M2"]^(3/2), type = "o" )
  
  
  # M4
  plot( x = n * Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 / sum((z - mean(z))^4), 
        y = EEPN[ ,4] ) # !
  
  A <- cbind( x = n * Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 / sum( (z - mean(z))^4 ), 
              y = EEPN[ ,4] )   # !!!!!!!
  A[ ,1] / A[ ,2]
  A
  
  cbind( Moments[ ,"M4_boot"], 
         EEPN[ ,4] * (sum( (z - mean(z))^4 ) / (n * Moments[ ,"M2"]^2)) )
  
  
  
  
  plot( x = sum( (z - mean(z))^4 ) * EEPN[ ,4], 
        y = Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2, type = "o"   )
  plot( x = (sum( (z - mean(z))^4 ) * EEPN[ ,4]) / Moments[ ,"M2"]^2, 
        y = Moments[ ,"M4_boot"], type = "o"   )
  plot( x = Moments[ ,"M4_boot"],
        y = EEPN[ ,4], type = "o" )
  
  
  plot( x = EEPN[ ,4] * sum( (z - mean(z))^4 ) / Moments[ ,"M2"],
        y = Moments[ ,"M4_boot"] )
  
  
  
  
  
  plot( x = Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 / (n * sum( (z - mean(z))^4 ) ), 
        y = EEPN[ ,4], type = "o" )
  
  
  
  a <- Moments[ ,"M4_boot"] * n
  b <- Moments[ ,"M2"]^2
  d <- sum( (z - mean(z))^4 )
  g <- EEPN[ ,4] - WBPN[ ,4]
  # g <- EEPN[ ,4]
  
  plot( x = a * b / d, y = g, type = "o" )
  plot( x = a, y = g / (b/d), type = "o" )
  plot( x = a, y = g * d / b, type = "o" )
  plot( x = a/n, y = (g * d / b) / n, type = "o" )
  plot( x = a/n, y = g * d / (n * b), type = "o", xlim = range(a/n), ylim = range(g * d / (n * b)) )
  
  head( cbind( a * b / d, y = g ), 70 )
  head( cbind( x = a/n, y = g * d / (n * b) ), 70 )
  cbind( x = a/n, y = g * d / (n * b) )
  cbind( x = a/n, y = g * d / (n * b) )[100, ]
  
  
  idx <- 10
  (EEPN[idx, 4] - WBPN[idx, 4]) * sum( (z - mean(z))^4 ) / (n * Moments[idx, "M2"]^2 )
  Moments[idx, "M4_boot"]
  
  
  M4BBFUN( x = z )
  M4BBFUN( x = z, M2 = Moments[idx, "M2"] )
  
  
  
  
  
  
  A <- cbind( sum( (z - mean(z))^4 ) * EEPN[ ,4], 
              Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2 )
  plot( x = A[ ,1], y = A[ ,2], type = "o" )
  head(A)
  tail(A)
  plot( A[ ,1] / A[ ,2] )
  plot( A[ ,1] - A[ ,2] )
  
  f <- A[ ,2] / A[ ,1]
  plot( A[ ,1] * f, y = A[ ,2], type = "o" )
  
  A[100, ]; f[100]
  A[100, 1] * f[100]; A[100, 2]
  WBPN[100, 4]
  Moments[100, "M4_boot"]
  sum( (z - mean(z))^4 ) / n * EEPN[100, 4]
  kurtFUN(z) * EEPN[100, 4]
  
  
  
  
  tmp <- cbind( sum( (z - mean(z))^4 ) * EEPN[ ,4], 
                Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2,
                sum( (z - mean(z))^4 ) * EEPN[ ,4] * f,
                sum( (z - mean(z))^4 ) * EEPN[ ,4] * f / Moments[ ,"M2"]^2,
                Moments[ ,"M4_boot"] )
  head(tmp)
  
  
  cbind( kurtFUN(x = z) * ( 1 + (EEPN[ ,4] - 1) / (1 - WBPN[ ,4]) ),
         Moments[ ,"M4_boot"] )
  
  
  
  
  
  w_bar_mat <- t( do.call( cbind, lapply( scl_vec, FUN = function(scl) { alpha * scl } ) ) )
  headleft(w_bar_mat)
  
  
  
  plot( x = 10^4 * Moments[ ,"M6_boot"] * Moments[ ,"M2"]^3 / sum((z - mean(z))^6), y = EEPN[ ,6], type = "o" ) # !
  
  cbind( Moments[ ,"M6_boot"] * Moments[ ,"M2"]^3,
         sum((z - mean(z))^6) * EEPN[ ,6] )
  
  
  cbind( Moments[ ,"M4_boot"] * Moments[ ,"M2"]^2,
         sum((z - mean(z))^4) * EEPN[ ,4] )
  
  
  
  plot( x = Moments[ ,"M6_boot"], y = sum((z - mean(z))^6) * EEPN[ ,6] / Moments[ ,"M2"]^3, type = "o" )
  plot( x = Moments[ ,"M6_boot"], y = (sum((z - mean(z))^6) * EEPN[ ,6]) / Moments[ ,"M2"]^3, type = "o" )
  plot( x = n * Moments[ ,"M6_boot"] * Moments[ ,"M2"]^3 / sum( (z - mean(z))^6 ), 
        y = EEPN[ ,6] )
  
  
  
  A <- cbind( x = (EEPN[ ,6]^(1/6))^5.66, 
              y = Moments[ ,"M6_boot"] * Moments[ ,"M2"]^3 )
  slope <- (A[1, "y"] - A[nrow(A), "y"]) / (A[1, "x"] - A[nrow(A), "x"])
  intercept <- A[nrow(A), "y"] - slope * A[nrow(A), "x"]
  plot( x = A[ ,"x"], y = A[ ,"y"] )
  abline(a = intercept, b = slope)
  
  
  
  
  
  
  
  
  
  
  
  