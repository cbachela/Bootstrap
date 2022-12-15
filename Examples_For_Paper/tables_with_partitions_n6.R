  
  
  ############################################################################
  ### EXACT MOMENTS - TABLES WITH PARTITIONS, N = 6
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.06.2022
  # First version:    16.06.2022
  # --------------------------------------------------------------------------
  
  
 
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # n out of n, n = 6
  # --------------------------------------------------------------------------
  
  n <- 6
 
  # Partitions
  
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # Permutations
  
  lP <- list()
  for( j in 1:ncol(parts_mat) ) {
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    lP[[j]] <- perm_unique
  }
  gridnodes_all <- do.call( rbind, lP )
  dim(gridnodes_all)
  sum(draws[2, ])
  comb( n = 2*n-1, k = n )
  
  
  
  
  round( cbind( t(parts_mat), t(draws), prob, prob_tot ), 5 )
  
  
  
  
  # --------------------------------------------------------------------------
  # m out of n, n = 6, m = 4
  # --------------------------------------------------------------------------

  n <- 6
  m <- 4
  
  parts_mat <- partitions::parts( m )
  parts_mat <- rbind( parts_mat, 
                      matrix(0, nrow = n-m, ncol = ncol(parts_mat)) )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  lP <- list()
  for( j in 1:ncol(parts_mat) ) {
    if ( length(unique(parts_mat[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat[ ,j])
    }
    lP[[j]] <- perm_unique
  }
  gridnodes_all <- do.call( rbind, lP )
  
  dim(gridnodes_all)
  sum(draws[2, ])
  # comb( n = 2*n-1, k = n )
  # comb( n = 2*m-1, k = m )
  # comb( n = 2*n-1, k = m )
  # comb( n = 2*m-1, k = n )
  # # ?????????
  
  round( cbind( t(parts_mat), t(draws), prob, prob_tot ), 5 )
  
  
    
  
  
  
  
  
  
  
  