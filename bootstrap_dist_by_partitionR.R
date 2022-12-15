  
  
  ############################################################################
  ### BOOTSTRAP - LP NORM AND ENTROPY BY PARTITION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     30.04.2021
  # First version:    16.04.2021
  # --------------------------------------------------------------------------
  
  
  # Compute distance measures among partition members
  # Find patterns between partition stats
  # Pairwise distances among group (partition) members
  # Check which partitions have the same distance measures
  # Expected distances
  # Compare with BB
  # Compare with block bootstrap
  # Cauchy-Schwarz
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(boot)
  require(volesti)
  require(RP)
  require(slolz)
  
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm_n80.rds") )
  n <- env$n
  B <- env$B
  gridnodes <- env$gridnodes
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  prob_dmult <- env$prob_dmult
  prob_tot <- env$prob_tot
  id_match <- env$id_match
  lp <- env$lp
  lp_relative <- env$lp_relative
  llp_pairwise <- env$llp_pairwise
  etpy <- env$etpy
  etpy_relative <- env$etpy_relative
  etpy_pairwise <- env$etpy_pairwise
  etpy_relative_ls <- env$etpy_relative_ls
  p_k <- env$p_k
  p_l <- env$p_l
  p_tot <- env$p_tot
  draws_emp <- env$draws_emp
  
  
  
  plot( head(prob_tot, 1000) )
  n
  class(B)
  dim(gridnodes)
  
  
  
  # Unique grid nodes
  gridnodes_unique <- gridnodes[!duplicated(gridnodes), ]
  
  # Level sets
  gridnodes_unique_sort <- t( apply( gridnodes_unique, 1, function(x) { rev(sort(x)) } ) )
  gridnodes_ls <- gridnodes_unique_sort[!duplicated(gridnodes_unique_sort), ]

  
  
  
  
  # --------------------------------------------------------------------------
  # Compute distance measures among partition members
  # --------------------------------------------------------------------------
  
  x_bar <- rep(1/n, n)
  id_match_unique <- unique(id_match)
  
  
  lGNU <- list()
  lGNUS <- list()
  for ( l in seq(along = id_match_unique) ) {
    
    idx <- which( id_match == id_match_unique[l] )
    lGNU[[l]] <- gridnodes_unique[idx, ]
    lGNUS[[l]] <- gridnodes_unique_sort[idx, ]
    
  }
  
  # Entropy
  etpyFUN <- function(P) 
  {
    if ( is.null(dim(P)) ) {
      P <- matrix( P, nrow = 1 )
    }
    apply( P, 1, entropy, exponential = FALSE, base = 2, eps = 1e-15 )
  }
  etpy <- lapply( lGNU, FUN = etpyFUN )
  etpy_s <- lapply( lGNUS, FUN = etpyFUN )
  
  
  # Relative entropy
  etpyrelFUN <- function(P) 
  {
    if ( is.null(dim(P)) ) {
      P <- matrix( P, nrow = 1 )
    }
    apply( P, 1, relentropy, q = x_bar, exponential = FALSE, base = 2, eps = 1e-15 )
  }
  etpy_rel <- lapply( lGNU, FUN = etpyrelFUN )
  etpy_rel_s <- lapply( lGNUS, FUN = etpyrelFUN )
  
  
  # Lp
  lpFUN <- function(P) 
  {
    if ( is.null(dim(P)) ) {
      P <- matrix( P, nrow = 1, ncol = length(P) )
    }
    p_vec <- c(0.5, 0.99, 1, 2, 5, 10)
    tmp <- lapply( p_vec,
                   FUN = function(x) { apply( P, 1, lpNorm, p = x ) } )
    ans <- do.call( cbind, tmp )
    colnames(ans) <- paste0("p=",p_vec)
    return( ans )
  }
  lp <- lapply( lGNU, FUN = lpFUN )
  lp_s <- lapply( lGNUS, FUN = lpFUN )
  
  
  # Relative Lp
  lprelFUN <- function(P) 
  {
    if ( is.null(dim(P)) ) {
      P <- matrix( P, nrow = 1, ncol = length(P) )
    }
    p_vec <- c(0.5, 0.99, 1, 2, 5, 10)
    tmp <- lapply( p_vec,
                   FUN = function(x) { apply( P, 1, lpNorm, 
                                              x_bar = x_bar, p = x ) } )
    ans <- do.call( cbind, tmp )
    colnames(ans) <- paste0("p=",p_vec)
    return( ans )
  }
  lp_rel <- lapply( lGNU, FUN = lprelFUN )
  lp_rel_s <- lapply( lGNUS, FUN = lprelFUN )
  
  
  etpy
  etpy_s
  etpy_rel
  etpy_rel_s
  lp
  lp_s
  lp_rel
  lp_rel_s
  
  
  n_etpy <- unlist( lapply( etpy, FUN = function(x) { length(unique(x)) } ) )
  n_etpy_s <- unlist( lapply( etpy_s, FUN = function(x) { length(unique(x)) } ) )
  n_etpy_rel <- unlist( lapply( etpy_rel, FUN = function(x) { length(unique(x)) } ) )
  n_etpy_rel_s <- unlist( lapply( etpy_rel_s, FUN = function(x) { length(unique(x)) } ) )
  
  n_lp <- t( do.call( cbind, lapply( lp, FUN = function(X) { apply( X, 2, function(x) { length(unique(x)) } ) } ) ) )
  n_lp_s <- t( do.call( cbind, lapply( lp_s, FUN = function(X) { apply( X, 2, function(x) { length(unique(x)) } ) } ) ) )
  n_lp_rel <- t( do.call( cbind, lapply( lp_rel, FUN = function(X) { apply( X, 2, function(x) { length(unique(x)) } ) } ) ) )
  n_lp_rel_s <- t( do.call( cbind, lapply( lp_rel_s, FUN = function(X) { apply( X, 2, function(x) { length(unique(x)) } ) } ) ) )
  
  unlist( lapply( lGNU, function(x) { ifelse( is.null(nrow(x)), 1, nrow(x)) } ) ) 
  # debugonce( lpFUN )
  lpFUN( P = lGNU[[ length(lGNU) ]] )
  
  range( n_etpy )
  range( n_etpy_s )
  range( n_etpy_rel )
  range( n_etpy_rel_s )
  
  range( n_lp )
  range( n_lp_s )
  range( n_lp_rel )
  range( n_lp_rel_s )
  
  
  etpy[[ which( n_etpy_rel > 1 ) ]]
  
  
 

  
  # --------------------------------------------------------------------------
  # Check which partitions have the same distance measures
  # --------------------------------------------------------------------------
  
  # Question:
  # Can we deterimine bounds on the fraction of distinct distance measures?
  
  
  # On analytic partitions
  
  etpy_parts <- apply( parts_mat / n, 2, entropy, exponential = TRUE, base = exp(1), eps = 1e-15 )
  etpy_rel_parts <- apply( parts_mat / n, 2, relentropy, q = x_bar, exponential = TRUE, base = exp(1), eps = 1e-15 )
  
  p_vec <- c(0.5, 0.99, 1, 2, 5, 10)
  lp_parts <- lapply( p_vec, FUN = function(p) { 
                      apply( parts_mat / n, 2, lpNorm, p = p ) } )
  lp_parts <- do.call( cbind, lp_parts )
  lp_rel_parts <- lapply( p_vec, FUN = function(p) { 
                          apply( parts_mat / n, 2, lpNorm, x_bar = x_bar, p = p ) } )
  lp_rel_parts <- do.call( cbind, lp_rel_parts )
  
  
  
  
  etpy_parts_unique <- unique(etpy_parts)
  length(etpy_parts_unique); ncol(parts_mat)
  length(etpy_parts_unique) / ncol(parts_mat)
  
  etpy_rel_parts_unique <- unique(etpy_rel_parts)
  length(etpy_rel_parts_unique); ncol(parts_mat)
  length(etpy_rel_parts_unique) / ncol(parts_mat)
  
  lp_parts_unique <- apply( lp_parts, 2, unique )
  lapply(lp_parts_unique, length); ncol(parts_mat)
  unlist(lapply(lp_parts_unique, length)) / ncol(parts_mat)
  
  lp_rel_parts_unique <- apply( lp_rel_parts, 2, unique )
  lapply(lp_rel_parts_unique, length); ncol(parts_mat)
  unlist(lapply(lp_rel_parts_unique, length)) / ncol(parts_mat)
  
  
  # -------------------
  # Expected distances
  # -------------------
  
  etpy_parts_mu <- sum(etpy_parts * prob_tot)
  etpy_rel_parts_mu <- sum(etpy_rel_parts * prob_tot)
  lp_parts_mu <- t(prob_tot) %*% lp_parts
  lp_rel_parts_mu <- t(prob_tot) %*% lp_rel_parts
  
  
  plot( density( etpy_parts ) ) 
  abline( v = etpy_parts_mu )
  abline( v = mean(unlist(etpy)), col = 2 )
  
  plot( density( etpy_rel_parts ) ) 
  abline( v = etpy_rel_parts_mu )
  abline( v = mean(unlist(etpy_rel)), col = 2 )
  
  plot( density( lp_parts[ ,4] ) ) 
  abline( v = lp_parts_mu[4] )
  abline( v = mean( unlist(lapply( lp, function(x) { x[ ,4] } )) ), col = 2 )
  
  
  
  
  N <- 10^4
  plot( x = etpy_parts[1:N], y = log(prob[1:N])  )
  plot( x = log(tail(prob, N)), y = tail(etpy_parts, N) )
  
  plot( x = log(prob_tot[1:N]), y = log(prob[1:N]) )
  plot( x = prob_tot[1:N], y = prob[1:N] )
  
  
  
  plot.ldensity( lapply( lp_parts_unique[-3], function(x) { density( (x - mean(x)) / sd(x) ) } ),
                 fillin = FALSE )   
  
  plot( density( lp_parts[ ,1] ) )
  
  plot( density( lp_parts[ ,4] ) )
  plot( density( lp_parts[ ,6] ) )
  

  
  # on empirical samples
  
  etpy_unique <- unique( unlist(etpy) )
  length(etpy_unique)
  length(etpy)  

  etpy_rel_unique <- unique( unlist(etpy_rel) )
  length(etpy_rel_unique)
  length(etpy_rel)  
  
    
  
  
  # --------------------------------------------------------------------------
  # Find patterns between partition stats
  # --------------------------------------------------------------------------
  
  N <- 10^3
  plot( x = log(draws[1, 1:N]), log(draws[2, 1:N]) )
  
  plot( x = log(draws[1, 1:N]), log(prob[1:N]) )
  plot( x = log(draws[2, 1:N]), log(prob_tot[1:N]) )
  plot( x = log(draws[1, 1:N]), etpy_parts[1:N] )
  plot( x = log(prob_tot[1:N]), etpy_parts[1:N] )
  plot( x = log(prob[1:N]), etpy_parts[1:N] )
  
  
  
  plot( log( draws[2, 1:N] ) )
  draws[ ,1:10]
  
  
  
  
  # --------------------------------------------------------------------------
  # Pairwise distances among group (partition) members
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # Compare with BB
  # --------------------------------------------------------------------------
  
  B <- nrow(gridnodes)
  S <- simplex( d = n )
  RPS <- .rpS( object = S, n_sim = B )
  wmat <- getSamples(RPS)
  
  etpy_wmat <- apply( wmat, 1, entropy, exponential = FALSE, base = 2, eps = 1e-15 )
  etpy_rel_wmat <- apply( wmat, 1, relentropy, q = x_bar, exponential = FALSE, base = 2, eps = 1e-15 )
  lp_wmat <- lapply( p_vec, function(p) { apply( wmat, 1, lpNorm, p = p ) } )  
  lp_wmat <- do.call( cbind, lp_wmat )
  lp_rel_wmat <- lapply( p_vec, function(p) { apply( wmat, 1, lpNorm, x_bar = x_bar, p = p ) } )  
  lp_rel_wmat <- do.call( cbind, lp_rel_wmat )

  etpy_gn <- apply( gridnodes, 1, entropy, exponential = FALSE, base = 2, eps = 1e-15 )
  etpy_rel_gn <- apply( gridnodes, 1, relentropy, q = x_bar, exponential = FALSE, base = 2, eps = 1e-15 )
  lp_gn <- lapply( p_vec, function(p) { apply( gridnodes, 1, lpNorm, p = p ) } )  
  lp_gn <- do.call( cbind, lp_gn )
  lp_rel_gn <- lapply( p_vec, function(p) { apply( gridnodes, 1, lpNorm, x_bar = x_bar, p = p ) } )  
  lp_rel_gn <- do.call( cbind, lp_rel_gn )
  
  
  plot( x = sort(etpy_gn), y = sort(etpy_wmat) )
  
  

  plot.ldensity( apply( cbind( etpy_gn, etpy_wmat ), 2, density ) )
  abline( v = etpy_parts_mu )
  plot.ldensity( apply( cbind( etpy_rel_gn, etpy_rel_wmat ), 2, density ) )
  abline( v = etpy_rel_parts_mu )
  
  i = 4
  plot.ldensity( apply( cbind( lp_gn[ ,i], lp_wmat[ ,i] ), 2, density ) )
  abline( v = lp_parts_mu[i] )
  
  i = 4
  plot.ldensity( apply( cbind( lp_rel_gn[ ,i], lp_rel_wmat[ ,i] ), 2, density ) )
  abline( v = lp_rel_parts_mu[i] )
  
  
  # Try to estimate dispersion of expected distances
  j <- 1
  d_parts <- lp_parts[ ,j]
  d_parts_mu <- lp_parts_mu[j]
  d_emp <- cbind( lp_gn[ ,j], lp_wmat[ ,j] )
  # d_emp <- cbind( lp_gn[ ,j], lp_gn[ ,j] )
  
    
  # P <- rdirichlet( n = 10^2, alpha = prob_tot )
  # mu <- P %*% d_parts
  
  # Since total probability strongly concentrates, it suffices to evaluate just
  # the most probable partitions
  idx <- which( cumsum(prob_tot) > 0.999 )[1]
  P2 <- rdirichlet( n = 10^2, alpha = prob_tot[1:idx] )
  mu2 <- P2 %*% d_parts[1:idx]  
  
  plot.ldensity( apply( d_emp, 2, density ) )
  abline( v = d_parts_mu )
  lines( density( mu ), col = 1 )
  lines( density( mu2 ), col = 2 )
  
  
  
  
  
  boxplot( etpy_parts )
  abline( h = mean(etpy_gn) )
  abline( h = mean(etpy_wmat), col = 2 )
  
  
  MCMCpack::ddirichlet( x = wmat[1, ], alpha = rep(1/n, n) )
    
  
  

    
  # --------------------------------------------------------------------------
  # Compare with block bootstrap
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # Cauchy-Schwarz
  # --------------------------------------------------------------------------
  
  set.seed(1234)
  x <- rnorm(n)
  z <- x^2
  
  l_x <- lpNorm(x = x, p = 2)
  Boot <- boot( data = x,
                statistic = meanBoot,
                R = 10^3 )  
  range(Boot$t)
  range(abs(Boot$t))
  range(l_x * lp_parts[ ,4])
  range(lp_parts[ ,4])  
  
  head(lp_parts)
  range( lp_parts[1:10^4 ,4] )
  
  theta_hat <- apply( gridnodes, 1, function(p) { sum(p * x) } )
  l_x * range( lp_gn[ ,4] )
  range(theta_hat)
  l_x * range( lp_gn[ ,4] )
  
  sum(x^2) * range(lp_gn[ ,4])^2
  range(theta_hat)^2  
    
  
  
  n <- 10^2
  B <- 10^5 + 1
  x <- rnorm(n)
  z <- x^2
  
  l_x <- lpNorm(x = x, p = 2)
  Boot <- boot( data = x,
                statistic = meanBoot,
                R = B )  
    
  S <- simplex( d = n )
  RPS <- .rpS( object = S, n_sim = B )
  wmat <- getSamples(RPS)
  
  lpw <- apply( wmat, 1, lpNorm, p = 2 )
  l_x * range(lpw)
  range(Boot$t)  
  range(x)  
  
  hw <- apply( wmat, 1, entropy, exponential = FALSE, base = 2, eps = 1e-15 )
  range(hw)
  
  
  range( wmat %*% x )  
  l_x / range(hw)
    
  
  
    
  
  
  
  
  
  