  
  
  ############################################################################
  ### BOOTSTRAP GRID NODES - Lp DISTANCES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.04.2021
  # First version:    16.04.2021
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(Rmpfr)
  require(partitions)
  require(boot)
  require(volesti)
  require(RP)
  require(slolz)
  
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  # Analytic results
  # Bootstrap sampling
  # Lp-Norms
  # Entropy
  # Pairwise Lp-distances
  # Pairwise relative entropy
  
  

  # --------------------------------------------------------------------------
  # Analytic results
  # --------------------------------------------------------------------------
  
  n <- 30
  B <- 10^3 + 1
  comb( n = 2*n-1, k = n )
  combMpfr( n = 2*n-1, k = n )
  combMpfr2( n = 2*n-1, k = n  )
  P(n)
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = TRUE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  
  # Reorder results by decreasing total probability
  ordering <- rev( order( prob_tot) ) 
  parts_mat <- parts_mat[ ,ordering]
  draws <- draws[ ,ordering]
  prob <- prob[ordering]
  prob_tot <- prob_tot[ordering]
  prob_dmult <- prob_dmult[ordering]
  
  
  plot( prob_tot[ 1:10^3 ] )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bootstrap sampling
  # --------------------------------------------------------------------------
  
  
  # Sample grid nodes
  tic <- Sys.time()
  gridnodes <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample(1:n, replace = TRUE)
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )
  
  # Unique grid nodes
  gridnodes_unique <- gridnodes[!duplicated(gridnodes), ]
  
  # Level sets
  gridnodes_unique_sort <- t( apply( gridnodes_unique, 1, function(x) { rev(sort(x)) } ) )
  gridnodes_ls <- gridnodes_unique_sort[!duplicated(gridnodes_unique_sort), ]
  
  
  dim( gridnodes )
  dim( gridnodes_unique )
  dim( gridnodes_ls )
  
  
  
  # --------------------------------------------------------------------------
  # Node activation (log) probabilities
  # --------------------------------------------------------------------------
  
  p_k <- apply( gridnodes_unique * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  p_k_sort <- apply( gridnodes_unique_sort * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  p_l <- apply( gridnodes_ls * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  
  draws_emp <- gridDraws( parts_mat = t(gridnodes_ls) * n )
  colnames(draws_emp) <- paste0("part", 1:ncol(draws_emp))
  p_tot <- p_l * draws_emp[2, ]
  
  sum( p_tot )
  
  plot( p_k - p_k_sort )  # same same
 
  
  
  # --------------------------------------------------------------------------
  # Match sampled gridnodes to analytic level sets, aka partitions
  # --------------------------------------------------------------------------
  
  id_match <- apply( gridnodes_unique_sort * n, 1, prodlim::row.match, 
                     table = t(as.matrix(parts_mat[ ,1:10^4])) )
  id_match_unique <- unique( id_match )

  any( is.na(id_match) )
  
  sum(which( id_match > 5000 )) / length(id_match)
  
  
  
  # What is the rank in total probability of sampled nodes
  prob_log <- log(prob)
  p_k_rnk <- lapply( rev(sort(log(p_k))), FUN = function(p) { which( prob_log == p )[1] } )

  prob_tot_log <- log( prob_tot ) 
  p_tot_rnk <- lapply( log(p_tot), FUN = function(p) { which( prob_tot_log == p )[1] } )
  
  unlist(p_tot_rnk)
  length(unlist(p_tot_rnk))
  head(unlist(p_tot_rnk))
  
  head(id_match)
  sum( id_match == unlist( p_tot_rnk ) )
  tmp <- which( id_match != unlist( p_tot_rnk ) )
  
  
  id_match[tmp]
  cbind( gridnodes_ls[tmp[1], ] * n, parts_mat[ ,id_match[tmp[1]]])
  
  
  p_tot[tmp[1]]
  prob_tot[id_match[tmp[1]]]
  
  
  cbind( p_tot, prob_tot[id_match] )
  plot( p_tot - prob_tot[id_match] )
  plot( p_tot - prob_tot[unlist(p_tot_rnk)] )
  plot( p_tot,  prob_tot[unlist(p_tot_rnk)] )
  plot( p_l,  prob[unlist(p_l_rnk)] )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Lp-Norms
  # --------------------------------------------------------------------------
  
  p_vec <- c(0.5, 1, 2, 10)
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
  
  headleft(lp)
  headleft(lp_relative)
  
  
  
  # --------------------------------------------------------------------------
  # Entropy
  # --------------------------------------------------------------------------
  
  x_bar <- rep(1/n, n)
  etpy <- apply( gridnodes_unique, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative <- apply( gridnodes_unique, 1, relentropy, q = x_bar, 
                          base = 2, exponential = FALSE )
  etpy_ls <- apply( gridnodes_ls, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative_ls <- apply( gridnodes_ls, 1, relentropy, q = x_bar, 
                             base = 2, exponential = FALSE )
  
  
  plot( etpy )
  plot( etpy_relative )
  
  
  
  # --------------------------------------------------------------------------
  # Pairwise Lp-distances
  # --------------------------------------------------------------------------
  
  pairwiseLPDist <- function( i, p = 1 ) 
  { 
    apply( gridnodes_unique[-i, ], 1, function(x) { (sum(abs(x - gridnodes_unique[i, ])^p))^(1/p) } )
  }
  
  llp_pairwise <- list()
  for ( i in seq(along = p_vec) ) {
    llp <- lapply( 1:nrow(gridnodes_unique), FUN = pairwiseLPDist, p = p_vec[i] )
    llp_pairwise[[i]] <- do.call( cbind, llp ) 
  }
  
  headleft(llp_pairwise[[1]] )
  
  
  
  # --------------------------------------------------------------------------
  # Pairwise relative entropy
  # --------------------------------------------------------------------------
  
  pairwiseEntropy <- function( i ) 
  { 
    FUN <- function(p) 
    { 
      relentropy( p = p, q = gridnodes_unique[i, ], base = 2, exponential = FALSE )
    }
    ans <- apply( gridnodes_unique[-i, ], 1, FUN = FUN )
    return( ans )
  }
  
  letpy <- lapply( 1:nrow(gridnodes_unique), FUN = pairwiseEntropy )
  etpy_pairwise <- do.call( cbind, letpy )
  
  headleft( etpy_pairwise )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$n <- n
  env$B <- B
  env$gridnodes <- gridnodes
  env$parts_mat <- parts_mat
  env$draws <- draws
  env$prob <- prob
  env$prob_dmult <- prob_dmult
  env$prob_tot <- prob_tot
  env$id_match <- id_match
  env$p_k <- p_k
  env$p_l <- p_l
  env$p_tot <- p_tot
  env$draws_emp <- draws_emp
  env$lp <- lp
  env$lp_relative <- lp_relative
  env$lp_ls <- lp_ls
  env$lp_relative_ls <- lp_relative_ls
  env$llp_pairwise <- llp_pairwise
  env$etpy <- etpy
  env$etpy_relative <- etpy_relative
  env$etpy_pairwise <- etpy_pairwise
  env$etpy_ls <- etpy_ls
  env$etpy_relative_ls <- etpy_relative_ls
 
  saveRDS( env, file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm_n80.rds") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm_n80.rds") )
    n <- env$n
    B <- env$B
    gridnodes <- env$gridnodes
    parts_mat <- env$parts_mat
    draws <- env$draws
    prob <- env$prob
    prob_dmult <- env$prob_dmult
    prob_tot <- env$prob_tot
    lp <- env$lp
    lp_relative <- env$lp_relative
    llp_pairwise <- env$llp_pairwise
    etpy <- env$etpy
    etpy_relative <- env$etpy_relative
    etpy_pairwise <- env$etpy_pairwise
    p_k <- env$p_k
    p_l <- env$p_l
    p_tot <- env$p_tot
    draws_emp <- env$draws_emp
    
    ### 
    prob <- exp(prob_dmult)
    prob_tot <- as.numeric(prob * draws[2, ])
    sum(prob_tot)
    ###
    
    
    
    # Convergence of total probability 
    plot( cumsum(prob_tot)[1:10^3] )
    abline(h = 0.99)
    
    idx <- which( cumsum( prob_tot ) > 0.99 )[1]
    idx
    
    
    # What is the rank in total probability of sampled nodes
    prob_sort_log <- log( rev( sort( prob ) ) ) 
    p_l_rnk <- lapply( rev(sort(log(p_l))), FUN = function(p) { which( prob_sort_log == p ) } )
    p_l_rnk[[1]]
    
    plot( prob_sort_log[1:1000] )
    points( rev( sort( log( p_l ) ) ), col = 2 )
    
    prob_tot_sort_log <- log( rev( sort( prob_tot ) ) ) 
    p_tot_rnk <- lapply( rev(sort(log(p_tot))), FUN = function(p) { which( prob_tot_sort_log == p )[1] } )
    
    plot( x = 1:length(p_tot_rnk), y = unlist(p_tot_rnk) )
    
    
    head( unlist(p_tot_rnk) )
    
    
    
    # Identify high total probability level sets
    ordering <- rev( order( prob_tot ) )
    
    
    plot( cumsum(prob_tot[ordering[1:length(p_l)] ]) )
    points( cumsum(rev(sort(p_tot))), col = 2 )
    
    tmp <- cbind( analytic = rev(sort(prob_tot))[1:length(p_l)],
                  empirical = rev(sort(p_tot)) )
    head(tmp)
    
    
    tmp1 <- cbind( analytic = prob[ordering[1:length(p_l)]],
                   empirical = p_l[rev(order(p_tot))] )
    head(tmp1)
    
    
    tmp2 <- cbind( analytic = draws[2, ordering[1:length(p_l)]],
                   empirical = draws_emp[2, rev(order(p_tot))] )
    head(tmp2)
    
    head( tmp1 * tmp2 )
    
    
    
    
    plot( log(tmp1), xlab = "analytic", ylab = "empirical" )
    
    plot( log(tmp1[ ,1]), ylim = range(log(tmp1)) )
    points( log(tmp1[ ,2]), col = 2 )
    
    plot( tmp2[ ,1], ylim = range(tmp2) )
    points( tmp2[ ,2], col = 2 )
    
    
    
    
    
    o <- ordering[1:length(p_k)]
    for ( i in 1:o ) {
      which( p_k == prob[o[i]] )
    }
    
    p_k %in% prob[o]
    
    barplot( head(cbind(p_k, prob[o]), 100), beside = TRUE, col = 1:2 )
    
    

   
    
    plot( cumsum( p_tot ) )

    
    
    
    
    apply( lp, 2, function(x) { length(unique(x)) } )
    apply( lp_relative, 2, function(x) { length(unique(x)) } )
    length(p_l)
    
    length(unique(etpy))
    length(unique(etpy_relative))
    
    
    plot(  x = log(p_l), y = etpy_relative_ls )
    points(  x = unique(p_l), y = unique(etpy_relative_ls), col = 2 )
    
    
    
    
    
    
    # Notes:
    # The number of unique Lp-norms for p<1 and p>10 equals the number 
    # of level sets
    # The number of unique relative entropy levels (to the centroid) equals
    # the number of level sets
    
    # Do they coincide? I.e., can level sets be identified with relative entropy?
    
    
    plot( log(p_k, base = 2), etpy_relative )
    
    
    id_match[ which( etpy_relative == etpy_relative[2] ) ]
    
    
    
    
    
    
    
    
    par(mfrow = c(2, 2))
    for ( i in 1:ncol(lp) ) {
      # boxplot( lp[ ,i] ) 
      plot( density( lp[ ,i] ) )
    }
    
    par(mfrow = c(2, 2))
    for ( i in 1:ncol(lp_relative) ) {
      # boxplot( lp_relative[ ,i] )
      plot( density( lp_relative[ ,i] ) )
    }
    
    
    
    dim(etpy_pairwise)
    
    
    plot( etpy )
    plot( etpy_relative )  
    
    
    plot( x = lp[ ,1], y = etpy )
    plot( x = lp[ ,3], y = etpy )
    plot( x = lp[ ,4], y = etpy )
    
    plot( x = lp_relative[ ,3], y = etpy_relative )
    
  
    
    boxplot( t(llp_pairwise[[1]]), beside = TRUE )
    boxplot( t(llp_pairwise[[2]]), beside = TRUE )
    boxplot( t(llp_pairwise[[3]]), beside = TRUE )
    boxplot( t(llp_pairwise[[4]]), beside = TRUE )
    
  
          
    
  }
  
  
  
  
