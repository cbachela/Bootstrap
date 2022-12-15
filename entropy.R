    
    
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


  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  n <- 30
  k <- n
  B <- 10^4 + 1
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  prob_dmult_log <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = TRUE )
  
  tmp <- cbind( prob = prob, prob_dmult = prob_dmult, idem = (prob == prob_dmult) )
  head(tmp, 1000)
  sum(tmp[ ,3]); nrow(tmp)
  
  ###
  prob <- prob_dmult
  prob_tot <- prob * draws[2, ]
  ###
  
  
  
  comb( n = 2*n-1, k = n )
  P(n)
  
  
  
  # Sample
  gridnodes <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
  }
  
  gridnodes_unique <- gridnodes[!duplicated(gridnodes), ]
  
  tmp <- data.table( t( apply(gridnodes, 1, function(x) { rev(sort(x)) } ) ) )
  gridnodes_ls <- tmp[!duplicated(tmp), ]
  
  # tbl_ls <- t( apply( gridnodes_unique * n, 1, function(x) { rev(sort(x)) } ) )
  # p_k_idx <- apply( tbl_ls, 1, prodlim::row.match, table = t(as.matrix(parts_mat)) )
  # p_k <- prob[ p_k_idx ]
  # 
  # tbl_ls <- t( apply( gridnodes_ls * n, 1, function(x) { rev(sort(x)) } ) )
  # p_l_idx <- apply( tbl_ls, 1, prodlim::row.match, table = t(as.matrix(parts_mat)) )
  # p_l <- prob[ p_l_idx ]
  
  p_k <- apply( gridnodes_unique * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  p_l <- apply( gridnodes_ls * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  draws_emp <- gridDraws( parts_mat = t(gridnodes_ls * n) )
  colnames(draws_emp) <- paste0("partition", 1:ncol(draws_emp))
  p_total <- p_l * draws_emp[2, ]

  sum( p_total )
  
  
  # Convergence of total probability 
  plot( cumsum(rev(sort(prob_tot)))[1:length(p_total)] )
  points( cumsum(rev(sort(p_total))), col = 2 )
  abline(h = 0.99)
  
  delta <- cumsum(rev(sort(prob_tot)))[1:length(p_total)] - cumsum(rev(sort(p_total)))
  plot( delta )
  
  idx <- which( cumsum( rev( sort( prob_tot ) ) ) > 0.99 )[1]
  idx
  
  
  
  plot( draws[2, ] )
  dim( draws )
  
  max( draws[2, ] )
  
  
  p_total %in% prob_tot
  p_total %in% (prob_dmult * draws[2, ])    # same same
  
  length(unique(p_total))
  length(p_l)
  
  
  
  plot( p_total )
  plot( density( draws[2, ] ) )
  plot( draws[2, ] )
  
  
  dim(draws)
  dim(gridnodes)
  dim(gridnodes_unique)
  dim(gridnodes_ls)
  
  length( unique(p_k) )
  length( unique(p_l) )
  length( p_l )
  length(p_k)
  
  
  
  
  
  
  # Lp norm
  
  p_vec <- c(0.25, 0.5, 0.75, 1, 2, 3, 5, 10, 20)
  x_bar <- rep(1/n, n)
  lp <- matrix( NA, nrow = nrow(gridnodes_unique), ncol = length(p_vec) )
  lp_relative <- lp
  lp_ls <- matrix( NA, nrow = nrow(gridnodes_ls), ncol = length(p_vec) )
  lp_relative_ls <- lp_ls
  for ( i in seq(along = p_vec) ) {
    lp[ ,i] <- apply( gridnodes_unique, 1, lpNorm, p = p_vec[i] )
    lp_relative[ ,i] <- apply( gridnodes_unique, 1, lpNorm, x_bar = x_bar, p = p_vec[i] )
    lp_ls[ ,i] <- apply( gridnodes_ls, 1, lpNorm, p = p_vec[i] )
    lp_relative_ls[ ,i] <- apply( gridnodes_ls, 1, lpNorm, x_bar = x_bar, p = p_vec[i] )
  }
  
  # Entropy
  x_bar <- rep(1/n, n)
  etpy <- apply( gridnodes_unique, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative <- apply( gridnodes_unique, 1, relentropy, q = x_bar, 
                          base = 2, exponential = FALSE )
  etpy_ls <- apply( gridnodes_ls, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative_ls <- apply( gridnodes_ls, 1, relentropy, q = x_bar, 
                             base = 2, exponential = FALSE )
  
  
  lp_unique <- apply( lp, 2, function(x) { length(unique(x)) } )
  lp_relative_unique <- apply( lp_relative, 2, function(x) { length(unique(x)) } )
  lp_ls_unique <- apply( lp_ls, 2, function(x) { length(unique(x)) } )
  lp_relative_ls_unique <- apply( lp_relative_ls, 2, function(x) { length(unique(x)) } )
  lp_unique
  lp_relative_unique
  lp_ls_unique
  lp_relative_ls_unique 
  length(p_l)
  length(unique(etpy_relative))
  
  i <- 1
  plot( x = lp[ ,i], y = p_k )
  plot( x = lp_relative[ ,i], y = p_k )
  
  
  plot(  x = unique(p_k), y = unique(etpy_relative) )
  plot(  x = unique(lp[ ,8]), y = unique(etpy_relative) )
  
  plot( x = log(p_l), y = etpy_relative_ls )
  plot( x = log(p_l), y = lp_ls[ ,1] )
  
  
  
  
  barplot( log(p_k) )
  barplot( etpy_relative )
  
  
  idx <- which(p_k == p_k[1])
  plot(  x = p_k[idx], y = etpy_relative[idx] )
  
  
  
  
  
  
  # Typicality
  
  H <- entropy( p = rep(1/k, k), exponential = FALSE, base = 2 )
  eps <- 1e-03
  a <- 2^(n * (-H + eps))
  b <- 2^(n * (-H - eps))
  
  # boxplot( log(p_k) )
  # abline( h = log(a) )
  # abline( h = log(b), col = 2 )
  range(log(p_l))
  log(a); log(b)
  
  boxplot( log(p_l) )
  abline( h = log(a) )
  abline( h = log(b), col = 2 )
  
  H; range( -1/n * log(p_l, base = 2))
  
  
  H; range( -1/n * log(prob, base = 2) )
  
  plot( -1/n * log(prob, base = 2) )
  abline(h = H)  

  
  plot( x = prob_tot, y = -1/n * log(prob, base = 2) )
  abline(h = H)  
  
  
  i <- which( -1/n * log(prob, base = 2) == H )
  i
  parts_mat[ ,i]  
  
  
  #// H corresponds to avg log-prob of vertex nodes. Why ???
  
  
  
  
