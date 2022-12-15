  
  
  ############################################################################
  ### BAYESIAN BOOTSTRAP SAMPLES - DISTANCE TO ORDINARY BOOTSTRAP GRID
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     21.04.2021
  # First version:    21.04.2021
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
  
  
  
  
  # --------------------------------------------------------------------------
  # Parameters
  # --------------------------------------------------------------------------
  
  n <- 10^4 * 5
  B <- 10^3 + 1
  set.seed(1234)
  z <- rnorm(n)^2
  p_vec <- c(0.5, 0.99, 1, 2, 5, 10, 20)
  
  
  # --------------------------------------------------------------------------
  # Sample Bayesian bootstrap samples
  # --------------------------------------------------------------------------
  
  S <- simplex( d = n )
  RPS <- .rpS( object = S, n_sim = B )
  wmat <- getSamples(RPS)
  
  
  # --------------------------------------------------------------------------
  # Sample ordinary bootstrap samples
  # --------------------------------------------------------------------------
  
  gridnodes <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
  }
  gridnodes_unique <- gridnodes[!duplicated(gridnodes), ]
  
  
  dim(gridnodes)
  dim(gridnodes_unique)
  
  
  
  # --------------------------------------------------------------------------
  # Lp distance from bayesian samples to gridnodes
  # --------------------------------------------------------------------------
  
  pairwiseLPDist <- function( i, p = 1 ) 
  { 
    apply( wmat[-i, ], 1, function(x) { (sum(abs(x - gridnodes_unique[i, ])^p))^(1/p) } )
  }
  
  llp_pairwise <- list()
  for ( i in seq(along = p_vec) ) {
    llp <- lapply( 1:nrow(wmat), FUN = pairwiseLPDist, p = p_vec[i] )
    llp_pairwise[[i]] <- do.call( cbind, llp ) 
  }
  
  
  
  boxplot( t(llp_pairwise[[1]]), beside = TRUE )
  boxplot( t(llp_pairwise[[2]]), beside = TRUE )
  boxplot( t(llp_pairwise[[3]]), beside = TRUE )
  boxplot( t(llp_pairwise[[4]]), beside = TRUE )
  
  
 
  
  # --------------------------------------------------------------------------
  # Entropy
  # --------------------------------------------------------------------------
  
  x_bar <- rep(1/n, n)
  etpy_w <- apply( wmat, 1, entropy, base = 2, exponential = FALSE )
  etpy_gn <- apply( gridnodes, 1, entropy, base = 2, exponential = FALSE )
  etpy_rel_w <- apply( wmat, 1, relentropy, 
                       q = x_bar, base = 2, exponential = FALSE )
  etpy_rel_gn <- apply( gridnodes, 1, relentropy, 
                        q = x_bar, base = 2, exponential = FALSE )
  
  
  plot( x = sort(etpy_rel_w), y = sort(etpy_rel_gn) )
  
  ldens <- apply( cbind( etpy_rel_w, etpy_rel_gn), 2, density )
  plot.ldensity( ldens )
  
  ldens <- apply( cbind( etpy_w, etpy_gn), 2, density )
  plot.ldensity( ldens )
  
  
  
  
  # --------------------------------------------------------------------------
  # LP
  # --------------------------------------------------------------------------
  
  p <- 10
  lp_w <- apply( wmat, 1, lpNorm, p = p )
  lp_gn <- apply( gridnodes, 1, lpNorm, p = p )
  lp_rel_w <- apply( wmat, 1, lpNorm, x_bar = rep(1/n, n), p = p )
  lp_rel_gn <- apply( gridnodes, 1, lpNorm, x_bar = rep(1/n, n), p = p )
  
  
  plot( x = sort(lp_w), y = sort(lp_gn) )
  
  ldens <- apply( cbind( lp_w, lp_gn), 2, density )
  plot.ldensity( ldens )
  
  
  ldens <- apply( cbind( lp_rel_w, lp_rel_gn), 2, density )
  plot.ldensity( ldens )
  
  lpNorm( x = rep(1/n, n), p = p )
  
  
  
  
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
 
  saveRDS( env, file = paste0(wd, "waRehouse/xxx.rds") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    env <- readRDS( file = paste0(wd, "waRehouse/xxx.rds") )
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
    
  
    
    
    
    ##############################################
    
    
    env_grid <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm.rds") )
    lp_grid <- env_grid$lp
    lp_grid_relative <- env_grid$lp_relative
    etpy_grid <- env_grid$etpy
    etpy_grid_relative <- env_grid$etpy_relative
    
    
    env_bb <- readRDS( file = paste0(wd, "waRehouse/bayesian_bootstrap_lp_norm.rds") )
    lp_bb <- env_bb$lp
    lp_bb_relative <- env_bb$lp_relative
    etpy_bb <- env_bb$etpy
    etpy_bb_relative <- env_bb$etpy_relative
    
    
    
    par(mfrow = c(2, 2))
    for ( i in 1:ncol(lp_grid) ) {
      ldens <- apply( cbind( lp_grid[ ,i], lp_bb[ ,i] ), 2, density )
      plot.ldensity( ldens, fillin = FALSE )
    }
    
    par(mfrow = c(2, 2))
    for ( i in 1:ncol(lp_grid_relative) ) {
      ldens <- apply( cbind( lp_grid_relative[ ,i], lp_bb_relative[ ,i] ), 2, density )
      plot.ldensity( ldens, fillin = FALSE )
    }
    
    
    dev.off()
    
    ldens <- apply( cbind( etpy_grid, etpy_bb ), 2, density )
    plot.ldensity( ldens, fillin = FALSE )
    
    ldens <- apply( cbind( etpy_grid_relative, etpy_bb_relative ), 2, density )
    plot.ldensity( ldens, fillin = FALSE )
    
    
          
    
  }
  
  
  
  
