  
  
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
  
  wd <- "H:/Papers Cyril/PhD_Papers/Bootstrap/"
  source( paste0(wd, "R/custom_functions.R") )
  
  
  # Research question:
  # Is it true that for n large enough, bootstrap weights belong to the 
  # typical set?
  # How can we show that a sampled point is typical?
  # What we know:
  # - the total probability of the typical set is close to one.
  # --> hence, p_k * draws[2, ] should be close to one.
  
  
  
  
  
  
  # Notes:
  # Probability level sets do not correspond to entropy level sets (some of the 
  # level sets have the same entropy). 
  # Also, entropy levels of probability level set members are not unique at all.
  
  
  
  
  
  
  
  # Differential entropy of Dirichlet symmetric random variable 
  diffentropyDirichlet <- function( alpha )
  {
    if ( all(alpha == 1) ) {
      H <- log( 1 / gamma( n ) )
    } else {
      n <- length(alpha)
      a0 <- sum(alpha)
      B <- prod( gamma(alpha) ) / gamma( a0 )
      H <- log(B) + (a0 - n) * digamma(a0) - sum( (alpha - 1) * digamma(alpha) )
    }
    return( H )
  }
  
  
  
  
  n <- 20
  B <- 10^4
  set.seed(1234)
  z <- rnorm(n)^2
  
  
  # Partitions
  tic <- Sys.time()
  parts_mat <- partitions::parts(n)
  (toc_parts <- Sys.time() - tic)
  dim(parts_mat)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_tot <- prob * draws[2, ]
  
  sum( prob * draws[2, ] )
  
  
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
  
  
  
  plot( p_k, H_LS )
  
  
  length(lWmat); length(H_unique); P(n); dim(parts_mat)[2]
  
  dim(DT)
  dim(DT_UGN)
  dim(DT_LS)
  length(H_unique)
  length(H_LS)
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
  
  
  # Typical set denoted A_{\epsilon}^{(n)} w.r.t. p(x) is the set of sequences 
  # x^n \in \mathcal{X}^n that satisfy the following:
  # 2^{???n(H(x)+\epsilon)} \leq p(x_1, . . . , x_n) \leq 2^{???n(H(x)???\epsilon))}
  # Or, equivalently:
  # A_{\epsilon}^{(n)} = \{ x_n \in X^n | p(x^n) \in [2^{-n(H(x)+\epsilon)}, 2^{???n(H(x)???\epsilon)}]
  
  
  eps <- 0.01
  2^(-n * (H0 + eps))
  2^(-n * (H0 - eps))
  range(p_k)
  
  
  # AEP
  range( -1/n * log(p_k, base = 2) )   # --> to H(X) in probability
  H0
  H1
  H2
  
  
  
  plot( x = prob, y = H_parts_mat, pch = 19, cex = 1.5, ylab = "entropy" )
  points( x = p_k, y = H_LS, col = 2, pch = 19 ) # , xlim = c(min(p_k), quantile(p_k, 0.01)) )
  # points( x = prob, y = H_LS[na.omit(match(H_parts_mat, H_LS))], col = 3, pch = 19 )
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
  env$prob <- prob
  env$DT_UGN <- DT_UGN     # unique grid nodes
  env$DT_LS <- DT_LS       # level sets
  env$H_UGN <- H_UGN
  env$H_LS <- H_LS
  env$H_unique <- H_unique # unique entropy levels
  env$p_k_idx <- p_k_idx
  env$p_k <- p_k           # level set probabilities
  
  saveRDS( env, file = paste0(wd, "waRehouse/bootstrap_grid_typical_set_n20.rds") )
  
  

  
  
  
