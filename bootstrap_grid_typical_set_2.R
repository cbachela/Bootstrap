  
  
  ############################################################################
  ### BOOTSTRAP GRID NODES - ENTROPY, TYPICAL SET AND AEP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.04.2021
  # First version:    16.04.2021
  # --------------------------------------------------------------------------
  
  
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
  
  
  
  # Research question:
  # Is it true that for n large enough, bootstrap weights belong to the 
  # typical set?
  
  
  
  # SAMPLING PROBABILITIES
  # ENTROPY
  # XXX
  # SAVE
  # ANALYZE
  

  
  
  # --------------------------------------------------------------------------
  # SAMPLING PROBABILITIES
  # --------------------------------------------------------------------------
  
  n <- 50
  B <- 10^7
  set.seed(1234)
  z <- rnorm(n)^2
  
  
  # Number of gridpoints
  tic <- Sys.time()
  parts_mat <- partitions::parts(n)
  (toc_parts <- Sys.time() - tic)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_tot <- prob * draws[2, ]
  
  
  comb( n = 2*n-1, k = n )
  P(n)
  exp(pi * sqrt(2*n/3)) / (4*n*sqrt(3)) # Harding Ramanujan
  
  
  # Bootstrap
  tic <- Sys.time()
  wmat <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )
  
  
  # Generate table of permutations (level sets)
  lwmat <- list()
  for ( i in 1:nrow(wmat) ) {
    lwmat[[i]] <- wmat[i, ]
  }
  tbl_list <- lwmat[ !duplicated(lapply(lwmat, sort)) ]
  level_sets <- matrix( unlist(tbl_list), ncol = n, byrow = TRUE )
  
  # Generate table of unique gridpoints
  lwmat <- list()
  for ( i in 1:nrow(wmat) ) {
    lwmat[[i]] <- wmat[i, ]
  }
  lGridpoints <- lwmat[ !duplicated( lwmat ) ]
  gp_unique <- do.call( rbind, lGridpoints )
  
  # Generate vector of probabilities associeted to each level set (p^(k))
  lP <- lapply( lGridpoints, FUN = function(x) { apply(parts_mat, 2, function(p) { all(p - rev(sort(x * n)) == 0) } ) } )
  p_k <- unlist( lapply( lP, FUN = function(p) { prob[ p ] } ) )
  p_k
  
  
  # Check if there are 2d-1 \choose d possible gridpoints
  comb( n = 2*n-1, k = n ); nrow(gp_unique)
  P(n); nrow(level_sets)
  
  
  
  
  head(wmat * n)
  t(parts_mat)
  gridpoints * n
  draws
  
  
  
  
  
  
  
  theta_gp <- unlist( lapply( lGridpoints, FUN = function(x) { t(x) %*% z } ) )
  
  plot( x = theta_values_gp, y = p_k )
  
  
  # Compute distribution function
  theta_gp_sort <- sort(theta_gp)
  d_gp <- p_k * NA
  for ( i in seq(theta_gp) ) {
    idx <- which(theta_gp <= theta_gp_sort[i])
    d_gp[i] <- sum(p_k[idx])
  }
  
  plot( x = theta_gp_sort, y = d_gp )
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ENTROPY
  # --------------------------------------------------------------------------
  
  
  
  
  
  a <- 1
  H1 <- digamma( n * a + 1 ) - digamma( a + 1) # discrete symmetric case
  H <- diffentropyDirichlet( alpha = rep(1, n) )
  H1
  H
  H0 <- digamma( a ) - digamma( n )
  
  basis <- 2
  E_unique_gp <- apply( gp_unique, 1, function(p) { entropy(p = p, exponential = FALSE, base = basis) } )
  E_level_set <- apply( level_sets, 1, function(p) { entropy(p = p, exponential = FALSE, base = basis) } )
  E_parts_mat <- apply( t(apply(parts_mat, 2, function(x) { x / sum(x)})), 1, function(p) { entropy(p = p, exponential = FALSE, base = basis) } )
  
  E <- H_unique_gp
  range(E)
  
  
  eps <- 0.1
  2^(n * (H0 + eps))
  2^(n * (H0 - eps))
  range(p_k)
  
  
  # AEP
  range( -1/n * log(p_k, base = basis) )   # --> to H(X)
  H
  H1
  H0
  
  
  
  plot( x = prob, y = E_parts_mat, pch = 19, cex = 1.5 )
  points( x = p_k, y = E_unique_gp, col = 2, pch = 19 ) # , xlim = c(min(p_k), quantile(p_k, 0.01)) )
  points( x = prob, y = E_level_set[match(E_parts_mat, E_level_set)], col = 3, pch = 19 )
  abline( h = -H0 )
  abline( h = mean(E1), col = 2 )
  abline( v = 2^(n * (H0 + eps)), col = "lightblue" )
  abline( v = 2^(n * (H0 - eps)), col = "darkblue" )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  
  n <- 10^3
  n_sim <- 10^3
  a <- 1
  P <- rdirichlet( n = n_sim, alpha = rep(a, n) )
  dim(P)
  
  E <- apply( P, 1, function(p) { entropy( p, exponential = FALSE, base = 2 ) } )
  
  
  H0 <- digamma( a ) - digamma( n )
  H1 <- digamma( n * a + 1 ) - digamma( a + 1) # discrete symmetric case
  range(E)
  H0
  H1
  
  
  
  
  # --------------------------------------------------------------------------
  # XXX
  # --------------------------------------------------------------------------
  
  
  
  
  
  # Notes
  # Probability level sets do not correspond to entropy level sets (there
  # are fewer unique entropy levels than unique probability levels)
  
  
  
  n <- 50
  B <- 10^3
  set.seed(1234)
  z <- rnorm(n)^2
  
  
  # Partitions
  tic <- Sys.time()
  parts_mat <- partitions::parts(n)
  (toc_parts <- Sys.time() - tic)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_scl <- prob / sum(prob)
  
  
  comb( n = 2*n-1, k = n )
  P(n)
  exp(pi * sqrt(2*n/3)) / (4*n*sqrt(3)) # Harding Ramanujan
  
  
  
  
  # Bootstrap loop
  wmat <- matrix( 0, nrow = B, ncol = n )
  tic <- Sys.time()
  for ( i in 1:B ) {
    
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
    
  }
  ( toc_boot <- Sys.time() - tic )
  
  # Create data.table and count how many times 
  # each unique grid node has been sampled
  DT <- data.table(wmat)
  counts <- DT[ ,.N, by = names(DT)][ ,"N"]
  
  # Create table with unique grid nodes
  DT_UGN <- DT[!duplicated(DT), ]
  
  # Create table with level sets
  tmp <- data.table( t( apply(DT, 1, sort) ) )
  DT_LS <- tmp[!duplicated(tmp), ]
  
  # Compute entropy
  basis <- exp(1)
  H_UGN <- apply( DT_UGN, 1, entropy, exponential = FALSE, base = basis )
  H_LS <- apply( DT_LS, 1, entropy, exponential = FALSE, base = basis )
  H_parts_mat <- apply( t(apply(parts_mat, 2, function(x) { x / sum(x)})), 1, function(p) { entropy(p = p, exponential = FALSE, base = basis) } )
  H_unique <- unique(H_LS)
  
  # Match probabilities from the partitions to each sampled level set (p^(k))
  tbl_ls <- t( apply( DT_LS * n, 1, function(x) { rev(sort(x)) } ) )
  p_k_idx <- apply( tbl_ls, 1, prodlim::row.match, table = t(as.matrix(parts_mat)) )
  p_k <- prob[ p_k_idx ]
  
  sum( prob * draws[2, ] )
  sum( p_k * draws[2, p_k_idx] )
  
  
  # # Group gridnodes by probability level sets (which does not correspond to entropy level sets)
  # lWmat <- list()
  # for ( i in 1:length(p_k) ) {
  #   lWmat[[i]] <- DT_UGN[which( H == H_unique[i] ), ]
  # }
  # 
  # # Check if all group members have same entropy
  # lH <- lapply( lWmat, FUN = function(W) { apply(w, 1, entropy, exponential = FALSE, base = basis ) } )
  
  
  
  
  
  length(lWmat); length(H_unique); P(n); dim(parts_mat)[2]
  
  dim(DT)
  dim(DT_UGN)
  dim(DT_LS)
  length(H_LS)
  length(H_unique)
  length(H_UGN)
  
  
  
  
  
  
  
  # Information theoretic concepts
  a <- 1
  H0 <- digamma( a ) - digamma( n )
  H1 <- digamma( n * a + 1 ) - digamma( a + 1) # discrete symmetric case
  H2 <- diffentropyDirichlet( alpha = rep(1, n) )

 
  range(H_unique)
  H0
  H1
  H2
  
  H <- H0
  eps <- 0.1
  2^(n * (H + eps))
  2^(n * (H - eps))
  range(p_k)
  
  boxplot( p_k )
  
  
  
  # AEP
  dim(t(parts_mat))
  dim(DT_LS)
  sum( prob * draws[2, ] )
  sum( p_k * draws[2, p_k_idx] )  # --> to 1
  
  range( -1/n * log(p_k, base = 2) )   # --> to H(X)
  H0
  H1
  H2
  
  
  
  plot( x = prob, y = H_parts_mat, pch = 19, cex = 1.5, ylab = "entropy" )
  points( x = p_k, y = H_LS, col = 2, pch = 19 ) # , xlim = c(min(p_k), quantile(p_k, 0.01)) )
  # points( x = prob, y = H_LS[match(H_parts_mat, H_LS)], col = 3, pch = 19 )
  abline( h = -H0 )
  abline( h = mean(E), col = 2 )
  abline( v = 2^(n * (H0 + eps)), col = "lightblue" )
  abline( v = 2^(n * (H0 - eps)), col = "darkblue" )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$n <- n
  env$B <- B
  env$z <- z
  env$parts_mat <- parts_mat
  env$draws <- draws
  env$prob <- prob
  env$DT_UGN <- DT_UGN     # unique grid nodes
  env$DT_LS <- DT_LS       # level sets
  env$H_UGN <- H_UGN
  env$H_LS <- H_LS
  env$H_unique <- H_unique # unique entropy levels
  env$p_k <- p_k           # level set probabilities
  env$p_k_idx <- p_k_idx
  
  saveRDS( env, file = paste0(wd, "waRehouse/bootstrap_grid_typical_set_n50.rds") )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_typical_set_n50.rds") )
  n <- env$n
  B <- env$B
  z <- env$z
  DT_LS <- env$DT_LS
  DT_UGN <- env$DT_UGN
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  H_UGN <- env$H_UGN
  H_LS <- env$H_LS
  H_unique <- env$H_unique
  p_k <- env$p_k
  p_k_idx <- env$p_k_idx
  
  
  
  H <- entropy( p = prob, exponential = FALSE, base = 2 )
  H
  
  
  # Information theoretic concepts
  a <- 1
  H0 <- digamma( a ) - digamma( n )
  H1 <- digamma( n * a + 1 ) - digamma( a + 1) # discrete symmetric case
  H2 <- diffentropyDirichlet( alpha = rep(1, n) )
  
  H <- entropy( p = rep(1/n, n), exponential = FALSE, base = 2 )
  eps <- 1e-19
  2^(n * (H + eps))
  2^(n * (H - eps))
  range(p_k)
  
  
  # AEP
  dim(t(parts_mat))
  dim(DT_LS)
  sum( prob * draws[2, ] )
  sum( p_k * draws[2, p_k_idx] )  # --> to 1
  
  range( -1/n * log(p_k, base = 2) )   # --> to H(X)
  H0
  H1
  H2
  
  
  
  
  
  
  
  

  
  # --------------------------------------------------------------------------
  # Grid node activation probabilities follow a multinomial distribution.
  # --------------------------------------------------------------------------
  
  #// Note: also holds for bootstrap case where we sample k out of n, k < n.
  
  n <- 60
  k <- 6
  parts_mat <- partitions::parts(n)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_scl <- prob / sum(prob)

  p_multinom <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n) )   
  
  cbind( prob, p_multinom )  # same same
  range( prob - p_multinom )
  
  
  
  p <- rep(1/k, k)
  H <- -sum( p * log(p, base = 2) )
  H <- entropy( p = p, exponential = FALSE, base = 2 )
  H
  barplot( prob )  
  abline( h = 2^(n * (-H + 0)) )
  range(prob)
  
  
  eps <- 1e-03
  a <- 2^(n * (-H + eps))
  b <- 2^(n * (-H - eps))
  boxplot( log(prob) )
  abline( h = log(a) )
  abline( h = log(b), col = 2 )
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  n <- 80
  k <- 6
  B <- 10^3
  parts_mat <- partitions::parts(n)
  prob <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n) )  
  
  # Sample
  gridnodes <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
  }
  gridnodes_unique <- gridnodes[!duplicated(gridnodes), ]
  
  tmp <- data.table( t( apply(gridnodes, 1, sort) ) )
  DT_LS <- tmp[!duplicated(tmp), ]
  tbl_ls <- t( apply( gridnodes_unique * n, 1, function(x) { rev(sort(x)) } ) )
  p_k_idx <- apply( tbl_ls, 1, prodlim::row.match, table = t(as.matrix(parts_mat)) )
  p_k <- prob[ p_k_idx ]
  
  
  eps <- 1e-03
  a <- 2^(n * (-H + eps))
  b <- 2^(n * (-H - eps))
  boxplot( log(prob) )
  abline( h = log(a) )
  abline( h = log(b), col = 2 )
  
  
  
  
  
