  
  
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
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
 
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 9
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
  w_bar <- rep(1/n, n)
  lNorm <- list()
  lNorm2wbar <- list()
  m_vec <- c(0.2, 0.5, 1:M)
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
  
  lp_mat_rel <- do.call( cbind, lapply( lNorm2wbar, FUN = function(x) { unlist(lapply(x, FUN = unique)) } ) )
  lp_mat_rel
  
  
  
  mean(z)
  mean(unlist(lTheta))  # same same
  sum( unlist(lNAP) * unlist(lTheta) ) # same same
  
  # M2
  M2 <- var(z) * ( sum(prob_tot * lp_mat[ ,4]^2) - 1/n )
  M2
  (sum( (z - mean(z))^2 ) / n) * (1 + (sum(prob_tot * lp_mat[ ,4]^2) - 1) / (1 - lpNorm(w_bar, p = 2)^2))
  sum( (z - mean(z))^2 ) / n^2
  lSigma <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^2 ) / length(theta) } )
  sum( prob_tot * unlist(lSigma) )
  plot( x = lp_mat[ ,4]^2, y = unlist(lSigma) )
  
  
  # M4
  lMom4 <- lapply( lTheta, function(theta) { sum( (theta - mean(z))^4 ) / length(theta) } )
  sum( prob_tot * unlist(lMom4) )
  M4BFUN(x = z, scl = FALSE)

  
  (sum( (z - mean(z))^4 ) / n) * sum(prob_tot * lp_mat[ ,6]^4)
  (sum( (z - mean(z))^4 ) / n) * (1 + (sum(prob_tot * lp_mat[ ,6]^4) - 1) / (1 - lpNorm(w_bar, p = 4)^4) )
  
  
  
  plot( x = (lp_mat[ ,6])^4, y = unlist(lMom4) )
  
  
  
  
  
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
  
  
  
  
  
  
  
  
  