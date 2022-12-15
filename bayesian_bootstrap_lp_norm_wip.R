  
  
  ############################################################################
  ### BAYESIAN BOOTSTRAP - Lp DISTANCES
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
  
  
  
  # Bootstrap sampling
  # Lp-Norms
  # Entropy
  # Pairwise Lp-distances
  # Pairwise relative entropy
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap sampling
  # --------------------------------------------------------------------------
  
  n <- 80
  B <- 10^3 + 1
  
  # Sample grid nodes
  S <- simplex( d = n )
  RPS <- .rpS( object = S, n_sim = B )
  wmat_unique <- getSamples(RPS)
  
  
  
  # --------------------------------------------------------------------------
  # Lp-Norms
  # --------------------------------------------------------------------------
  
  p_vec <- c(0.5, 1, 2, 10)
  x_bar <- rep(1/n, n)
  lp <- matrix( NA, nrow = nrow(wmat_unique), ncol = length(p_vec) )
  lp_relative <- lp
  for ( i in seq(along = p_vec) ) {
    lp[ ,i] <- apply( wmat_unique, 1, lpNorm, p = p_vec[i] )
    lp_relative[ ,i] <- apply( wmat_unique, 1, lpNorm, x_bar = x_bar, p = p_vec[i] )
  }
  
  headleft(lp)
  headleft(lp_relative)
  
  
  
  # --------------------------------------------------------------------------
  # Entropy
  # --------------------------------------------------------------------------
  
  x_bar <- rep(1/n, n)
  etpy <- apply( wmat_unique, 1, entropy, base = 2, exponential = FALSE )
  etpy_relative <- apply( wmat_unique, 1, relentropy, q = x_bar, 
                          base = 2, exponential = FALSE )
  
  plot( etpy )
  plot( etpy_relative )
  
  plot( density(etpy) )
  plot( density(etpy_relative) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Pairwise Lp-distances
  # --------------------------------------------------------------------------
  
  pairwiseLPDist <- function( i, p = 1 ) 
  { 
    apply( wmat_unique[-i, ], 1, function(x) { (sum(abs(x - wmat_unique[i, ])^p))^(1/p) } )
  }
  
  llp_pairwise <- list()
  for ( i in seq(along = p_vec) ) {
    llp <- lapply( 1:nrow(wmat_unique), FUN = pairwiseLPDist, p = p_vec[i] )
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
      relentropy( p = p, q = wmat_unique[i, ], base = 2, exponential = FALSE )
    }
    ans <- apply( wmat_unique[-i, ], 1, FUN = FUN )
    return( ans )
  }
  
  letpy <- lapply( 1:nrow(wmat_unique), FUN = pairwiseEntropy )
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
  
  saveRDS( env, file = paste0(wd, "waRehouse/bayesian_bootstrap_lp_norm.rds") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    env <- readRDS( file = paste0(wd, "waRehouse/bayesian_bootstrap_lp_norm.rds") )
    lp <- env$lp
    lp_relative <- env$lp_relative
    llp_pairwise <- env$llp_pairwise
    etpy <- env$etpy
    etpy_relative <- env$etpy_relative
    etpy_pairwise <- env$etpy_pairwise
    
    
    
    
    dim(lp)
    boxplot( scale( lp, TRUE, FALSE ), beside = TRUE )
    boxplot( scale( lp_relative, TRUE, FALSE ), beside = TRUE )
    
    
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
  
  
  
  
