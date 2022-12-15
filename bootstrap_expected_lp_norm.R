  
  
  ############################################################################
  ### EXPECTED LP-NORM OF GRID
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     03.02.2022
  # First version:    03.02.2022
  # --------------------------------------------------------------------------
  
  
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Parameters
  # --------------------------------------------------------------------------
  
  B <- 10^6
  M <- 6
  n_vec <- 3:20
  set.seed(1)
  z_vec <- rnorm( n = max(n_vec) )^2
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
  for ( i in seq(along = n_vec) ) {
    
    n <- n_vec[i]
    
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
    
    lp <- apply( parts_mat / n, 2, function(p) { sum(p^2) } )
    sum(lp * prob_tot)
    2 / n - 1 / n^2
    
    # p-norm of partitions
    tmp <- lapply( 1:M, function(m) { apply( parts_mat / n, 2, function(x) { sum(x^m) } ) } )
    lp_parts_mat <- do.call( cbind, tmp )
    lp_parts_mat
    
    # Compare with bootstrap sampler
    wmat <- matrix( 0, nrow = B, ncol = n )
    for ( b in 1:B ) {
      idx <- sample( 1:n, replace = TRUE )
      tbl <- table(idx)
      wmat[b, as.numeric(names(tbl))] <- tbl / n
    }
    tmp <- lapply( 1:M, function(m) { apply( wmat, 1, function(x) { sum(x^m) } ) } )
    lp_parts_mat_boot <- do.call( cbind, tmp )
    lp_parts_mat_boot
    
  }
  
  
  apply( lp_parts_mat, 2, function(lp) { sum(lp * prob_tot) } )
  apply( lp_parts_mat_boot, 2, mean )
  
  
  
  # M2
  n <- 20
  m2_ml <- sum( (z_vec[1:n] - mean(z_vec[1:n]))^2 ) / n
  scl2 <- n / (n - 1)
  m2_unb <- m2_ml * scl2
  var(z_vec[1:n]) - m2_unb
  mu_2 <- m2_ml / n
  
  (mu_2 + m2_unb * 1/n) / m2_unb 
  2 / n - 1 / n^2
  (n-1) / n^2 + 1/n
  apply( lp_parts_mat, 2, function(lp) { sum(lp * prob_tot) } )[2]
  
  
  
  
  # M3
  n <- 20
  m3_ml <- sum( (z_vec[1:n] - mean(z_vec[1:n]))^3 ) / n
  scl3 <- n^2 / ((n-1) * (n-2))
  m3_unb <- m3_ml * scl3
  mu_3 <- m3_ml / n^2
  a <- 3 * (2/n - 1/n^2) * 1/n
  b <- 2 * 1/n^2
  
  (mu_3 + m3_unb * (a - b)) / m3_unb
  mu_3 / m3_unb + a - b
  apply( lp_parts_mat, 2, function(lp) { sum(lp * prob_tot) } )[3]
  expected_3norm3 <- 1 / (n^2 * scl3) + a - b  # !!!!!!!!!!!!
  expected_3norm3
  ((n-1) * (n-2)) / n^4 + a - b   #!!!
  
  
  1 / n^2
  scl3 * ( 1 / (n^2 * scl3) )
  
  
  
  
  
    
  
  
  
  
  