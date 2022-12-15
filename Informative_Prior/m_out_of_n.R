  
  
  ############################################################################
  ### EXACT MOMENTS - INFORMATIVE PRIOR - m OUT OF n SAMPLING
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.06.2022
  # First version:    16.06.2022
  # --------------------------------------------------------------------------
  
  
  # Classical Bootstrap, m out of n sampling
  # Tables in the paper are generated in script tables_with_partitions_n6.R

  
  
  
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
  # Data
  # --------------------------------------------------------------------------
  
  n <- 8
  m <- 5
  n_boot <- 10
  B <- 10^3 + 1
  
  set.seed(1111)
  df <- 5
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  
  
  2/n - 1/n^2
  2/m - 1/m^2
  
  
  M2Exact( z = z, exp2norm2 = 2/n - 1/n^2 )
  M2Analytic( z = z )
  M2Exact( z = z, exp2norm2 = 2/m - 1/m^2 )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Run sampler (classical bootstrap), n out of n
  # --------------------------------------------------------------------------
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  boot_mat_grid <- boot_mat
  lP_boot <- list()
  lNorm_boot <- list()
  prob <- rep(1/n, n)
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample( 1:n, replace = TRUE, prob = prob )
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    }
    lP_boot[[j]] <- gridnodes
    lNorm_boot[[j]] <- apply( lP_boot[[j]], 1, lpNorm, p = 2 )
    boot_mat_grid[ ,j] <- gridnodes %*% z
  }
  
  
  mean( apply( boot_mat, 2, var ) )
  mean( apply( boot_mat_grid, 2, var ) )
  mean( unlist( lapply( lNorm_boot, function(x) { mean(x^2)} ) ) ); 2/n - 1/n^2
  
  
  # --------------------------------------------------------------------------
  # Run sampler (classical bootstrap), m out of n
  # --------------------------------------------------------------------------
  
  prob <- rep(1/n, n)
  boot_mat_grid_m <- matrix( NA, nrow = B, ncol = n_boot )
  lP_boot_m <- list()
  lNorm_boot_m <- list()
  for ( j in 1:n_boot ) {
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample( 1:n, size = m, replace = TRUE, prob = prob )
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / m
    }
    lP_boot_m[[j]] <- gridnodes
    lNorm_boot_m[[j]] <- apply( lP_boot_m[[j]], 1, lpNorm, p = 2 )
    boot_mat_grid_m[ ,j] <- gridnodes %*% z
  }
  
  
  mean( apply( boot_mat_grid_m, 2, var ) )
  mean( unlist( lapply( lNorm_boot_m, function(x) { mean(x^2)} ) ) ); 2/m - 1/m^2 # differnt
  
  
  
  
 
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact solution with partitions and permutations, n out of n
  # --------------------------------------------------------------------------
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # Permutations
  lP <- list()
  lTheta <- list()
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
    theta <- apply( perm_unique / n, 1, function(p) { sum(p * z) } )
    lTheta[[j]] <- theta
  }
  gridnodes_all <- do.call( rbind, lP )
  dim(gridnodes_all)
  sum(draws[2, ])
  comb( n = 2*n-1, k = n )
  
  # 2-norm of partitions
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  
  sum( prob_tot * lp_parts^2 )
  2/n - 1/n^2 # same same
  mean( unlist( lapply( lNorm_boot, function(x) { mean(x^2)} ) ) )
  
  
  
  
  # --------------------------------------------------------------------------
  # Exact solution with partitions and permutations, m out of n
  # --------------------------------------------------------------------------
  
  # Partitions
  parts_mat_m <- partitions::parts( m )
  parts_mat_m <- rbind( parts_mat_m, 
                      matrix(0, nrow = n-m, ncol = ncol(parts_mat_m)) )
  draws_m <- gridDraws( parts_mat = parts_mat_m, cnames = FALSE, mpfr = FALSE )
  prob_m <- draws2prob( draws = draws_m )
  prob_tot_m <- prob_m * draws_m[2, ]
  
  # Permutations
  lP_m <- list()
  lTheta_m <- list()
  for( j in 1:ncol(parts_mat_m) ) {
    if ( length(unique(parts_mat_m[ ,j])) > 1 ) {
      perm <- gtools::permutations( n = n, 
                                    r = n, 
                                    v = parts_mat_m[ ,j], 
                                    set = FALSE,
                                    repeats.allowed = FALSE )
      perm_unique <- perm[!duplicated(perm), ]
    } else {
      perm_unique <- t(parts_mat_m[ ,j])
    }
    lP_m[[j]] <- perm_unique
    theta <- apply( perm_unique / m, 1, function(p) { sum(p * z) } )
    lTheta_m[[j]] <- theta
  }
  gridnodes_all_m <- do.call( rbind, lP_m )
  dim(gridnodes_all_m)
  sum(draws_m[2, ])

  # 2-norm of partitions
  lp_parts_m <- apply( parts_mat_m / m, 2, lpNorm, p = 2 )
  
  sum( prob_tot_m * lp_parts_m^2 )
  mean( unlist( lapply( lNorm_boot_m, function(x) { mean(x^2)} ) ) )
  
  2/n - 1/n^2
  2/m - 1/m^2
  
  
  apply( gridnodes_all / n, 2, mean )
  apply( gridnodes_all_m / m, 2, mean )
  
  M2Exact( z = z, exp2norm2 = 2/n - 1/n^2 )
  M2Analytic( z = z, method = "classic" )
  
  mean( apply( boot_mat, 2, var ) )
  mean( apply( boot_mat_grid, 2, var ) ) 
  
  
  
  mean( apply( boot_mat_grid_m, 2, var ) ) 
  
  
  
  
  lLP2 <- list()
  lLP3 <- list()
  lLP4 <- list()
  lm_vec <- list()
  n_vec <- 1:50
  m_vec <- 1:max(n_vec)
  for ( i in seq(along = n_vec) ) {
    n_tmp <- n_vec[i]
    lm_vec[[i]] <- m_vec
    lp2 <- numeric(length(m_vec))
    lp3 <- lp2
    lp4 <- lp2
    for ( j in seq(along = m_vec) ) {
      m_tmp <- m_vec[j]
      parts_mat_tmp <- partitions::parts( m_tmp )
      if ( m_tmp < n_tmp ) {
        parts_mat_tmp <- rbind( parts_mat_tmp,
                                matrix(0, nrow = n_tmp - m_tmp, ncol = ncol(parts_mat_tmp)) )
      }
      draws_tmp <- gridDraws( parts_mat = parts_mat_tmp, cnames = FALSE, mpfr = FALSE )
      prob_tmp <- draws2prob( draws = draws_tmp )
      prob_tot_tmp <- prob_tmp * draws_tmp[2, ]
      lp_parts_tmp <- apply( parts_mat_tmp / m_tmp, 2, lpNorm, p = 2 )
      lp2[j] <- sum( prob_tot_tmp * lp_parts_tmp^2 )
      lp_parts_tmp <- apply( parts_mat_tmp / m_tmp, 2, lpNorm, p = 3 )
      lp3[j] <- sum( prob_tot_tmp * lp_parts_tmp^3 )
      lp_parts_tmp <- apply( parts_mat_tmp / m_tmp, 2, lpNorm, p = 4 )
      lp4[j] <- sum( prob_tot_tmp * lp_parts_tmp^4 )
    }
    lLP2[[i]] <- lp2
    lLP3[[i]] <- lp3
    lLP4[[i]] <- lp4
  }
  

  
  colors <- fBasics::divPalette(n = length(lm_vec), "PuOr")
  plot( x = lm_vec[[length(lm_vec)]], y = lLP2[[length(lLP2)]], ylim = c(0, max(lLP2[[length(lLP2)]])), type = "l" )
  for ( i in 1:length(lm_vec) ) {
    points( x = lm_vec[[i]], y = lLP2[[i]], col = colors[i], type = "o" )
  }
  
  plot( x = lm_vec[[length(lm_vec)]], y = lLP3[[length(lLP3)]], ylim = c(0, max(lLP3[[length(lLP3)]])), type = "l" )
  for ( i in 1:length(lm_vec) ) {
    points( x = lm_vec[[i]], y = lLP3[[i]], col = colors[i], type = "o" )
  }
  
  plot( x = lm_vec[[length(lm_vec)]], y = lLP4[[length(lLP4)]], ylim = c(0, max(lLP4[[length(lLP4)]])), type = "l" )
  for ( i in 1:length(lm_vec) ) {
    points( x = lm_vec[[i]], y = lLP4[[i]], col = colors[i], type = "o" )
  }

  
  
  m_vec <- unique( unlist( lm_vec) ) 
  mat2 <- matrix( NA, nrow = max(m_vec), ncol = length(lLP), 
                  dimnames = list(1:max(m_vec), 1:length(lLP)) )
  mat3 <- mat2
  mat4 <- mat2
  for ( i in 1:ncol(mat) ) {
    idx <- lm_vec[[i]]
    mat2[idx, i] <- lLP2[[i]]
    mat3[idx, i] <- lLP3[[i]]
    mat4[idx, i] <- lLP4[[i]]
  }  
  colnames(mat2) <- paste0("n=", 1:ncol(mat2) )
  rownames(mat2) <- paste0("m=", 1:nrow(mat2) )
  colnames(mat3) <- colnames(mat4) <- colnames(mat2)
  rownames(mat3) <- rownames(mat4) <- rownames(mat2)
  
  headleft(mat2, 11)
  headleft(mat3, 11)
  headleft(mat4, 11)
  
  
  
  
  plot( mat2[2, ] )
  plot( mat2[ ,1] )
  plot( mat3[2, ] )
  plot( mat3[ ,1] )
  

  
  
  FUN <- function(m, n) { (n + m - 1) / (n * m) }   #// = E_{P_{G_m}}(||\omega||_2^2)
  
  
  check_mat2 <- mat2 * NA
  mat2_tmp <- mat2 * NA
  check_mat2_tmp <- mat2 * NA
  for ( j in 1:ncol(mat2) ) {
    for ( i in 1:nrow(mat2) ) {
      m <- as.numeric(substr(rownames(mat2)[i], 3, nchar(rownames(mat2)[i])))
      n <- as.numeric(substr(colnames(mat2)[j], 3, nchar(colnames(mat2)[j])))
      check_mat2[i, j] <- FUN( m = m, n = n )
      if ( m <= n ) {
        check_mat2_tmp[i, j] <- FUN( m = m, n = n )  
        mat2_tmp[i, j] <- mat2[i, j]
      }
    }
  }

  delta2 <- mat2 - check_mat2
  delta2_tmp <- mat2_tmp - check_mat2_tmp
  
  plot( na.omit(as.numeric(delta2)) )
  plot( na.omit(as.numeric(delta2_tmp)) )
  
  
  
  
  # FUN <- function(m, n) { (5*m*n - 6*m + 2) / (m^2*n^2) }
  FUN <- function(m, n) { (n-1)*(n-2)/(n^2*m^2) + 3*(n+m-1)/(m*n^2) - 2/n^2 }
  
  check_mat3 <- mat3 * NA
  mat3_tmp <- mat3 * NA
  check_mat3_tmp <- mat3 * NA
  for ( j in 1:ncol(mat3) ) {
    for ( i in 1:nrow(mat3) ) {
      m <- as.numeric(substr(rownames(mat3)[i], 3, nchar(rownames(mat3)[i])))
      n <- as.numeric(substr(colnames(mat3)[j], 3, nchar(colnames(mat3)[j])))
      check_mat3[i, j] <- FUN( m = m, n = n )
      if ( m <= n ) {
        check_mat3_tmp[i, j] <- FUN( m = m, n = n )  
        mat3_tmp[i, j] <- mat3[i, j]
      }
    }
  }
  
  delta3 <- mat3 - check_mat3
  delta3_tmp <- mat3_tmp - check_mat3_tmp
  
  headleft(delta3)
  diag(delta3)
  
  

  # FUN <- function(n) { ((n-1)*(n-2))/n^4 + 4/n^2 - 3/n^3 }
  FUN <- function(n) { (5*n^2 - 6*n + 2) / n^4 }
  FUNBB <- function(n) { 6 / ((n+1)*(n+2)) }
  
  FUN(1:10)
  FUNBB(1:10)
  
  
  

  # write.csv( mat, file = "Q:/10_Mitarbeiterordner/CB/mat2.csv" )
  
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  