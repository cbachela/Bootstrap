  
  
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
  
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  

  # Bootstrap sampling
  # Lp-Norms
  # Entropy
  # Pairwise Lp-distances
  # Pairwise relative entropy
  
  

  # --------------------------------------------------------------------------
  # Bootstrap sampling
  # --------------------------------------------------------------------------
  
  n <- 20
  B <- 10^3 + 1
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  comb( n = 2*n-1, k = n )
  P(n)
  
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
  tmp <- t( apply( gridnodes_unique * n, 1, function(x) { rev(sort(x)) } ) )
  gridnodes_ls <- tmp[!duplicated(tmp), ]
  
  
  dim( gridnodes )
  dim( gridnodes_unique )
  dim( gridnodes_ls )
  
  
  
  # --------------------------------------------------------------------------
  # Node activation (log) probabilities
  # --------------------------------------------------------------------------
  
  p_k <- apply( gridnodes_unique * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  p_l <- apply( gridnodes_ls * n, 1, dmultinom, prob = rep(1/n, n), log = FALSE )
  
  draws <- gridDraws( parts_mat = t(gridnodes_ls) * n )
  colnames(draws) <- paste0("part", 1:ncol(draws))
  p_total <- p_l * draws[2, ]
  
  sum( p_total )
  headleft(draws)
 
  
  
  
  
  # --------------------------------------------------------------------------
  # Lp-Norms
  # --------------------------------------------------------------------------
  
  p_vec <- c(0.5, 1, 2, 10)
  x_bar <- rep(1/n, n)
  lp <- matrix( NA, nrow = nrow(gridnodes_unique), ncol = length(p_vec) )
  lp_relative <- lp
  for ( i in seq(along = p_vec) ) {
    lp[ ,i] <- apply( gridnodes_unique, 1, lpNorm, p = p_vec[i] )
    lp_relative[ ,i] <- apply( gridnodes_unique, 1, lpNorm, p = p_vec[i], x_bar = x_bar )
  }
  
  headleft(lp)
  headleft(lp_relative)

  x <- gridnodes_unique[1, ]
  x                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  lpNorm(x, p = 2)^2 - lpNorm(x_bar, p = 2)^2
  lpNorm( x = x, p = 2, x_bar = x_bar )^2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  
  
  c2 <- apply( parts_mat / n, 2, lpNorm, p = 2 )^2
  a2 <- lpNorm(x_bar, p = 2)^2
  b2 <- c2 - a2
  b2_check <- apply( parts_mat / n, 2, lpNorm, p = 2, x_bar = x_bar )^2
  cbind(b2, b2_check) # same same
  
  
  sum(prob_tot * c2) - 1/n
  sum(prob_tot * b2) # same same
  
  1 + (sum(prob_tot * c2) - 1) / (1 - 1/n)
  1/n  # same same
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Entropy
  # --------------------------------------------------------------------------
  
  x_bar <- rep(1/n, n)
  etpy <- apply( gridnodes_unique, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative <- apply( gridnodes_unique, 1, relentropy, q = x_bar, 
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
  env$lp <- lp
  env$lp_relative <- lp_relative
  env$llp_pairwise <- llp_pairwise
  env$etpy <- etpy
  env$etpy_relative <- etpy_relative
  env$etpy_pairwise <- etpy_pairwise
 
  saveRDS( env, file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm.rds") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm.rds") )
    n <- env$n
    B <- env$B
    lp <- env$lp
    lp_relative <- env$lp_relative
    llp_pairwise <- env$llp_pairwise
    etpy <- env$etpy
    etpy_relative <- env$etpy_relative
    etpy_pairwise <- env$etpy_pairwise
    
    
    
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
  
  
  
  
