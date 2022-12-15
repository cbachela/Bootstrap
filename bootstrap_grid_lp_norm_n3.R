  
  
  ############################################################################
  ### BOOTSTRAP GRID NODES - Lp DISTANCES BY LEVEL SET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.11.2021
  # First version:    17.11.2021
  # --------------------------------------------------------------------------
  
  
  
  # Gridnodes within probability level sets share they same 
  # - lp-norm for all p > 0
  # - lp-norm to centroid of simplex for all p > 0
  
  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(boot)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
 
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 5
  M <- 6
  n_boot <- 50
  B <- 10^4 + 1
  
  set.seed(1111)
  df <- 4
  x <- rt( n = n, df = df )
  z <- x^2

  
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
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
  comb( n = 2*n-1, k = n )
  
  
  # Weights for each theta
  lW <- lapply( 1:length(lP), FUN = function(i) { rep(prob[i], nrow(lP[[i]])) } )
  w <- unlist(lW)
  sum(w)
  
  # p-norm of grid nodes
  w_bar <- rep(1/n, n)
  lNorm <- list()
  lNorm2wbar <- list()
  m_vec <- 1:M
  for ( i in seq(along = m_vec) ) {
    m <- m_vec[i]
    lNorm[[i]] <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = m) } )
    lNorm2wbar[[i]] <- lapply( lP, FUN = function(P) { apply(P / n, 1, lpNorm, p = m, x_bar = w_bar) } )
  }
  
  
  # Node activation  probabilities
  param_multinomial <- rep(1/n, n)
  lNAP <- lapply( lP, FUN = function(P) { apply(P, 1, dmultinom, prob = param_multinomial, log = FALSE) } )
  lNAP
  nap <- unlist( lapply( lNAP, FUN = function(x) { x[1] } ) )
  cbind( unlist( lapply(lNAP, unique) ), prob )  # same same
  
  
  lp_mat <- do.call( cbind, lapply( lNorm, FUN = function(x) { unlist( lapply(x, FUN = unique) ) } ) )
  lp_mat
  
  # lp_mat_rel <- do.call( cbind, lapply( lNorm2wbar, FUN = function(x) { unlist(lapply(x, FUN = unique)) } ) )
  lp_mat_rel <- do.call( cbind, lapply( lNorm2wbar, FUN = function(x) { unlist(lapply(x, FUN = function(w) { w[1] } ) ) } ) )
  lp_mat_rel
  
  
  
  mean(z)
  mean(unlist(lTheta))  # same same
  sum( unlist(lNAP) * unlist(lTheta) ) # same same
  
  
  
  # Rimoldini_2013 Weighted skewness and kurtosis unbiased by sample size
  
  V1 <- sum(w_bar)
  V2 <- sum(w_bar^2)
  V3 <- sum(w_bar^3)
  V4 <- sum(w_bar^4)
  
  
  # M1
  lMu <- lapply( lTheta, FUN = function(theta) { mean(theta) } )
  sum( prob_tot * unlist(lMu) )
  mean(z)
  mean(z) * (1 + (sum(prob_tot * lp_mat[ ,3]^1) - 1) )
  # mean(z) * sum( prob_tot * lp_mat_rel[ ,3] )
  
  

  
  
  # M2
  M2 <- var(z) * ( sum(prob_tot * lp_mat[ ,2]^2) - 1/n )
  M2
  (sum( (z - mean(z))^2 ) / (n - 1)) * ( sum(prob_tot * lp_mat_rel[ ,2]^2) )
  (sum( (z - mean(z))^2 ) / n) * ( sum(prob_tot * lp_mat_rel[ ,2]^2) / (1 - lpNorm(w_bar, p = 2)^2) )
  (sum( (z - mean(z))^2 ) / n) * (1 + (sum(prob_tot * lp_mat[ ,2]^2) - 1) / (1 - lpNorm(w_bar, p = 2)^2))
  lSigma <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^2 ) / length(theta) } )
  sum( prob_tot * unlist(lSigma) )
  sum( (z - mean(z))^2 ) / n^2
  (sum( (z - mean(z))^2 ) / n) * lpNorm(w_bar, p = 2)^2  # !!!
  
  #// --> ||\bar{\omega}||_2^2 = 1/n
  
  
  
  scl2 <- V1^2 / (V1^2 - V2)
  scl2
  n / (n - 1)  # same same
  1 / (1 - lpNorm(w_bar, p = 2)^2)  # same same
  # w_tmp <- runif(n)
  # 1 / (1 - lpNorm(w_tmp / sum(w_tmp), p = 2)^2)  # not the same as with w_bar
  
  sum( (z - mean(z))^2 ) / n * scl2
  var(z)  
  
  sum( (z - mean(z))^2 ) / n^2
  sum( (z - mean(z))^2 ) / n * sum(prob_tot * lp_mat_rel[ ,2]^2) / (1 - lpNorm(w_bar, p = 2)^2)
  sum( (z - mean(z))^2 ) / n * scl2 * sum(prob_tot * lp_mat_rel[ ,2]^2)
  
  
  
  
  m2 <- sum( (z - mean(z))^2 ) / n
  a <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^2 ) } )
  A <- m2 * scl2 * a
  sum( unlist(lW) * A )
  M2
  
  idx <- cumsum( unlist(lapply(lW, length)) )
  # lIdx <- list()
  # for ( i in 1:(length(idx)-1) ) {
  #   if ( i == 1 ) {
  #     lIdx[[i]] <- 1:idx[i]
  #   } else {
  #     lIdx[[i]] <- (idx[i-1] + 1):idx[i]
  #   }
  # }
  # lIdx
  sum( prob_tot * A[idx] )
  
  
  # Check Appendix A
  a <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^2 ) } )
  sum( unlist(lW) * a )
  sum( prob_tot * lp_mat[ ,2]^2 ) - lpNorm(w_bar, p = 2)^2  # same same
  sum( prob_tot * (lp_mat[ ,2]^2 - lpNorm(w_bar, p = 2)^2) ) # same same
  
  sum( prob_tot * lp_mat[ ,2]^2 ); 2 / n - 1 / n^2  # same same
  
  
  
  
  
  
  
  # M3
  sum((z - mean(z))^3) / n
  lMom3 <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^3 ) / length(theta) } )
  sum( prob_tot * unlist(lMom3) )
  sum((z - mean(z))^3) / n^3
  sum((z - mean(z))^3) / n * lpNorm(w_bar, p = 3)^3  # !!!
  
  #// --> ||\bar{\omega}||_3^3 = 1/n^2
  
  
  # Check
  scl3 <- V1^3 / (V1^3 - 3*V1*V2 + 2*V3)
  scl3
  (n^2 / ((n-1) * (n-2)))  # same same
  
  sum( (z - mean(z))^3 ) / n * scl3
  sum( (z - mean(z))^3 ) / n * (n^2 / ((n-1) * (n-2)))  # same same
  
  
  M3BFUN( x = z, M2 = 1, scl = FALSE )
  sum( (z - mean(z))^3 ) / n^3
  sum( (z - mean(z))^3 ) / n * sum(prob_tot * lp_mat_rel[ ,3]^3) / (1 - lpNorm(w_bar, p = 3)^3)
  
  
  
  plot( x = lp_mat_rel[ ,5]^3, y = unlist(lMom3), type = "o" )
  plot( x = (sum( (z - mean(z))^3 ) * n / ((n-1) * (n-2))) * lp_mat_rel[ ,3]^3, 
        y = unlist(lMom3), type = "o" )
  
  
  
  
  ###
  w_bar <- rep(1/n, n)
  m3 <- sum( (z - mean(z))^3 ) / n
  # a <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^3 ) / (1 - 3*V2 + 2*V3) } )
  # a <- apply( gridnodes_all, 1, function(p) { sum( abs(p/n - w_bar)^3 ) } ) * scl3 # nope!
  a <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^3 ) } ) * scl3 # yes!
  A <- m3 * a
  M3BFUN( x = z, M2 = 1, scl = FALSE )
  sum( prob_tot * unlist(lMom3) )  # same same
  m3 * lpNorm(w_bar, p = 3)^3 # same same !!!
  sum( unlist( lW ) * A ) # same same
  ###
  
  idx <- cumsum( unlist(lapply(lW, length)) )
  sum( prob_tot * A[idx] )
  cbind(unlist(lMom3), A[idx])
  
  
  # Check Appendix A
  a1 <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^3 ) } )
  a2 <- apply( gridnodes_all, 1, function(p) { sum((p/n)^3) - 3*sum((p/n)^2) * sum(w_bar^2) + 2 * sum(w_bar^3) } )
  sum( unlist(lW) * a1 )
  sum( unlist(lW) * a2 )
  sum( prob_tot * lp_mat[ ,3]^3 ) - 3 * sum( prob_tot * lp_mat[ ,2]^2 ) * sum(w_bar^2) + 2 * sum(w_bar^3)
  
  
  
 
  
  
  
  a <- apply( parts_mat, 2, function(p) { sum( (p/n - w_bar)^3 ) } ) * scl3
  cbind( a / scl3, lp_mat_rel[ ,3]^3 )
  A <- m3 * a
  plot( A, unlist(lMom3), type = "o" )   # Works!
  tmp <- cbind( A, unlist(lMom3) )
  ordering <- order(tmp[ ,1])
  barplot( t(tmp[ordering, ]), beside = TRUE, col = 1:ncol(tmp) )
  tmp[ordering, ] * 100

  
  
  
  
  
  # M4
  lMom4 <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^4 ) / length(theta) } )
  m2 <- sum( (z - mean(z))^2 ) / n
  m4 <- sum( (z - mean(z))^4 ) / n
  denom <- ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
  scl4a <- (V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4))
  scl4b <- (3*V1^2 * (2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4))
  a <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^4 ) } )
  A <- (scl4a / denom * m4 - scl4b / denom * m2^2) * a
  
  a2 <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^2 ) } )
  a4 <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^4 ) } )
  Atmp <- scl4a * a4 / denom * m4 - scl4b * a4 / denom * m2^2

  
  idx <- cumsum( unlist(lapply(lW, length)) )
  sum( prob_tot * A[idx] )
  sum( prob_tot * Atmp[idx] )
  tmp <- cbind(unlist(lMom4), A[idx])
  # tmp <- cbind(unlist(lMom4), abs(Atmp[idx]) )
  barplot( t(tmp), beside = TRUE )
  cbind( tmp , tmp[ ,1] / tmp[ ,2] )
  
  plot( tmp[ ,1], tmp[ ,2] )
  
  
  M4BFUN(x = z, scl = FALSE)
  sum( prob_tot * unlist(lMom4) )  # same same
  sum( (z - mean(z))^4 ) / n^4 + 3 * ( sum( (z - mean(z))^2 )^2 ) * (1/n^4 - 1/n^5) # same same
  sum( unlist( lW ) * A ) # slightly different
  sum( unlist( lW ) * Atmp )
  
  
  # Check
  a1 <- apply( gridnodes_all, 1, function(p) { sum( (p/n - w_bar)^4 ) } )
  a2 <- apply( gridnodes_all, 1, function(p) { sum((p/n)^4) - 3*sum(w_bar^4) - 4*sum((p/n)^3)*sum(w_bar^2) + 6*sum((p/n)^2)*sum(w_bar^3) } ) 
  sum( unlist(lW) * a1 )
  sum( unlist(lW) * a2 )
  sum( prob_tot * lp_mat[ ,4]^4 ) - 3 * sum(w_bar^4) - 4*sum( prob_tot * lp_mat[ ,3]^3 ) * sum(w_bar^2) + 6*sum( prob_tot * lp_mat[ ,2]^2 ) * sum(w_bar^3)
  
  
  
  
  # Check: compare with unweighted
  scl4a_unw <- n * (n^2 - 2*n + 3) / ((n - 1)*(n - 2)*(n - 3))
  scl4b_unw <- 3*n * (2*n - 3) / ((n - 1)*(n - 2)*(n - 3))
  scl4a_unw * m4 - scl4b_unw * m2^2
  scl4a / denom * m4 - scl4b / denom * m2^2  # same same
  
  
  
  
  
  
  a <- apply( parts_mat, 2, function(p) { sum( (p/n - w_bar)^4 ) } )
  A <- (scl4a * m4 - scl4b * m2^2) * a
  plot( A, unlist(lMom4), type = "o" )  
  tmp <- cbind( A, unlist(lMom4) )
  ordering <- order(tmp[ ,1])
  barplot( t(tmp[ordering, ]), beside = TRUE, col = 1:ncol(tmp) )
  tmp[ordering, ] * 100
  
  plot( x = unlist(lMom4), y = a )
  plot( (tmp[ ,1] / tmp[ ,2]) )
  
  
  
  
  
  
  
  
  
  # M5
  lMom5 <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^5 ) / length(theta) } )
  lpNorm(w_bar, p = 5)^5
  1 / n^4
  #// --> ||\bar{\omega}||_5^5 = 1/n^4
  
  m2 <- sum( (z - mean(z))^2 ) / n
  m3 <- sum( (z - mean(z))^3 ) / n
  m5 <- sum( (z - mean(z))^5 ) / n
  10 * m3 * m2 / n^3 + (m5 - 10 * m3 * m2) / n^4
  sum( prob_tot * unlist(lMom5) )
  sum( (z - mean(z))^5 ) / n^5  # different (for n > 4) !!!
  
  
  
  
  
  # M7
  lMom7 <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^7 ) / length(theta) } )
  sum( prob_tot * unlist(lMom7) )
  sum( (z - mean(z))^7 ) / n^7
  (sum( (z - mean(z))^7 ) / n) * lpNorm(w_bar, p = 7)^7
  
  
  
  
  
  
  
  ##############
  lKurt <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^4 ) / (length(theta) * M2^2) } )
  sum( prob_tot * unlist(lKurt) )
  M4BFUN(x = z); M4BFUN(x = z, M2 = M2)
  
  
  
  A <- cbind( x = var(z) * (lp_mat[ ,4]^2 - 1/n), 
              y = unlist(lSigma) )
  plot( x = A[ ,"x"], y = A[ ,"y"] ) 
  var(z) * (lp_mat[ ,4]^2 - 1/n) - unlist(lSigma)
  
  
  A <- cbind( x = (lp_mat[ ,4]), 
              y = unlist(lSigma) * 10)
  plot( x = A[ ,"x"], y = A[ ,"y"], ylim = range(A[ ,"y"]), xlim = range(A[ ,"y"]) )
  points( x = A[ ,"y"], y = A[ ,"y"], col = 2 )
  points( x = A[ ,"x"]^2, y = A[ ,"y"], col = 3 )
  
  A[ ,"x"]^2 - A[ ,"y"]
  
  
  
  
  A <- cbind( x = (lp_mat[ ,6]) , 
              y = unlist(lKurt) / 400 )
  plot( x = A[ ,"x"], y = A[ ,"y"], ylim = range(A[ ,"y"]), xlim = range(A[ ,"y"]) )
  points( x = A[ ,"y"], y = A[ ,"y"], col = 2 )
  points( x = A[ ,"x"]^4, y = A[ ,"y"], col = 3 )
  A
  
  A[ ,"x"]^4 - A[ ,"y"]
  
  
  
  
  
  ks <- kurtFUN( x = z ) * (var(z) * (n-1) / n)^(4/2)
  ks / n + 3 - 3 / n
  sum( (z - mean(z))^4 ) / n^2 + 3 - 3 / n
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  
  tic <- Sys.time()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  mean( apply( boot_mat, 2, var ) )
  
 
  
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
  
  
  