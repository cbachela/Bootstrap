  
  
  ############################################################################
  ### BOOTSTRAP GRID NODES - VISUALIZE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.04.2021
  # First version:    16.04.2021
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(partitions)
  require(boot)
  require(volesti)
  require(RP)
  require(slolz)
  require(scales)
  
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  

  
  
  # --------------------------------------------------------------------------
  # n = 3
  # --------------------------------------------------------------------------
  
  n <- 3
  B <- 10^3
 
  # Number of gridpoints
  tic <- Sys.time()
  parts_mat <- partitions::parts(n)
  (toc_parts <- Sys.time() - tic)
  draws <- gridDraws(parts_mat, mpfr = FALSE)
  prob <- draws2prob(draws)
  prob_tot <- prob * draws[2, ]
  
  # Number of unique grid nodes
  comb( n = 2*n-1, k = n )
  
  # Number of unique node activation probabilities (level sets)
  P(n)
  
  # Bootstrap
  tic <- Sys.time()
  wmat <- matrix(0, nrow = B, ncol = n, dimnames = list(NULL, paste0("w", 1:n)))
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )
 
  # Create table with unique grid nodes
  UGN <- wmat[!duplicated(wmat), ]
  
  # Sample uniformly from the simplex
  S <- Simplex$new( d = n )
  samples <- S$runif( 10^5 * 3 )
  colnames(samples) <- colnames(wmat)
  
  
  # Plot
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], 
          xlab <- colnames(samples)[1], ylab = colnames(samples)[2], zlab = colnames(samples)[3],
          size = 1, col = "grey",
          box = FALSE, axes = TRUE, alpha = 0.4 )
  points3d( x = UGN[ ,1], y = UGN[ ,2], z = UGN[ ,3], size = 8 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # n = 4
  # --------------------------------------------------------------------------
  
  n <- 4
  B <- 10^4
  
  # Number of gridpoints
  tic <- Sys.time()
  parts_mat <- partitions::parts(n)
  (toc_parts <- Sys.time() - tic)
  draws <- gridDraws(parts_mat)
  prob <- draws2prob(draws)
  prob_scl <- prob / sum(prob)
  
  # Number of unique grid nodes
  comb( n = 2*n-1, k = n )
  
  # Number of unique node activation probabilities (level sets)
  P(n)
  
  # Bootstrap
  tic <- Sys.time()
  wmat <- matrix(0, nrow = B, ncol = n, dimnames = list(NULL, paste0("w", 1:n)))
  for ( i in 1:B ) {
    idx <- sample( 1:n, replace = TRUE )
    tbl <- table(idx)
    wmat[i, as.numeric(names(tbl))] <- tbl / n
  }
  ( toc_boot <- Sys.time() - tic )
  
  # Create table with unique grid nodes
  UGN <- wmat[!duplicated(wmat), ]
  
  # Sample uniformly from the simplex
  S <- Simplex$new( d = n )
  samples <- S$runif( 10^5 )
  samples <- fevBias(samples, q = 2)
  samples_fev <- fevBias(samples, q = 20)
  samples <- rbind(samples, samples_fev)
  colnames(samples) <- colnames(wmat)
  
  
  # Plot
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], 
          xlab <- colnames(samples)[1], ylab = colnames(samples)[2], zlab = colnames(samples)[3],
          size = 1, col = "grey", 
          box = FALSE, axes = TRUE, alpha = 0.5 )
  points3d( x = UGN[ ,1], y = UGN[ ,2], z = UGN[ ,3], size = 8 )
  
  
  
  
  
  
  
  
  
  
  