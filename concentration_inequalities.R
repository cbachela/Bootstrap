  
  
  ############################################################################
  ### BOOTSTRAP - CONCENTRATION INEQUALITIES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     24.04.2021
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
  
  
  
  
  n <- 6
  B <- 10^4 + 1
  set.seed(1111)
  x <- rnorm(n)
  z <- x^2
  comb( n = 2*n-1, k = n )
  P(n)
  
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
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
  
  
  
  # env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm_n80.rds") )
  # n <- env$n
  # B <- env$B
  # gridnodes <- env$gridnodes
  # parts_mat <- env$parts_mat
  # draws <- env$draws
  # prob <- env$prob
  # prob_dmult <- env$prob_dmult
  # prob_tot <- env$prob_tot
  # id_match <- env$id_match
  # lp <- env$lp
  # lp_relative <- env$lp_relative
  # llp_pairwise <- env$llp_pairwise
  # etpy <- env$etpy
  # etpy_relative <- env$etpy_relative
  # etpy_pairwise <- env$etpy_pairwise
  # etpy_relative_ls <- env$etpy_relative_ls
  # p_k <- env$p_k
  # p_l <- env$p_l
  # p_tot <- env$p_tot
  # draws_emp <- env$draws_emp
  
  
  
  
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
  dim(gridnodes_unique)

  
  # ########################
  # env_gn <- new.env()
  # env_gn$gridnodes_unique <- gridnodes_unique
  # saveRDS( env_gn, file = "H:/R/notingit/Bootstrap/waRehouse/gridnodes80.rds" )
  # env_gn <- readRDS( file = "H:/R/notingit/Bootstrap/waRehouse/gridnodes80.rds" )
  # ########################
  
  
  # Level sets
  gridnodes_unique_sort <- t( apply( gridnodes_unique, 1, function(x) { rev(sort(x)) } ) )
  gridnodes_ls <- gridnodes_unique_sort[!duplicated(gridnodes_unique_sort), ]
  
  nrow(gridnodes)
  nrow(gridnodes_unique)
  
  
  # --------------------------------------------------------------------------
  # Match sampled gridnodes to analytic level sets, aka partitions
  # --------------------------------------------------------------------------
  
  # id_match <- apply( gridnodes_unique_sort * n, 1, prodlim::row.match, 
  #                    table = t(as.matrix(parts_mat[ ,1:min(ncol(parts_mat), 10^4)])) )
  # id_match_unique <- sort( unique( id_match ) )
  # length(id_match_unique)
  # 
  # 
  # # Alternatively:
  # 
  # id_match <- apply( gridnodes_unique_sort * n, 1, 
  #                    function(x) { which( (colSums(parts_mat[ ,1:10^3] == x) == n) ) } )
  # id_na <- which( unlist( lapply( id_match, FUN = function(x) { identical(x, integer(0)) } ) ) )
  # if ( length(id_na) > 0 ) {
  #   id_match_na <- apply( gridnodes_unique_sort[id_na, ] * n, 1, 
  #                         function(x) { which( (colSums(parts_mat[ ,(10^3+1):10^5] == x) == n) ) } )
  #   id_match[ id_na ] <- id_match_na
  #   id_match <- unlist( id_match )
  # }
  
 
  
  hashFun <- function(x) { sum( (1/x^2)^3 ) }
  A <- parts_mat[ ,1:min(ncol(parts_mat), 10^4)]
  A[ A == 0 ] <- 0.5
  hash_parts <- apply( A, 2, hashFun )
  B <- gridnodes_unique_sort * n
  B[ B == 0 ] <- 0.5
  hash_gridnodes <- apply( B, 1, hashFun )
  
  id_match <- lapply( hash_gridnodes, FUN = function(x) { which( hash_parts == x ) } )
  id_na <- which( unlist( lapply( id_match, FUN = function(x) { identical(x, integer(0)) } ) ) )
  if ( length(id_na) > 0 ) {
    A2 <- parts_mat[ ,(10^4+1):10^5]
    A2[ A2 == 0 ] <- 0.5
    hash_parts_2 <- apply( A2, 2, hashFun )
    id_match_na <- lapply( hash_gridnodes[ id_na ], FUN = function(x) { which( hash_parts_2 == x ) } )
    id_match[ id_na ] <- unlist( id_match_na ) + 10^4
  }
  id_match <- unlist( id_match )
  id_match_unique <- sort( unique( id_match) )

  # Check
  check <- lapply( 1:length(id_match), 
                   FUN = function(i) { all( parts_mat[ ,id_match[i]] == gridnodes_unique_sort[i, ] * n ) } )
  range(check)
  which( check == 0 )
  
  
  
  
  # --------------------------------------------------------------------------
  # Compute distance measures on partitions
  # --------------------------------------------------------------------------
  
  
  ent_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, entropy, exponential = FALSE, base = exp(1) )
  lp_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 2 )
  lp_parts3 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 3 )
  lp_parts5 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 5 )
  lp_parts10 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 10 )
  lp_parts100 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 100 )
  lp_parts05 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 0.5 )
  
  
  
  
  # --------------------------------------------------------------------------
  # Sample from each partition
  # --------------------------------------------------------------------------
  
  N <- 10^5
  theta_mat <- matrix( NA, nrow = N, ncol = ncol(parts_mat) ) 
  P <- matrix( NA, nrow = N, ncol = nrow(parts_mat) )
  for ( k in 1:ncol(parts_mat) ) {
    for ( j in 1:N ) {
      P[j, ] <- sample( x = parts_mat[ ,k], replace = FALSE )
    }  
    theta_mat[ ,k] <- (P / n) %*% z
  }
  theta_mat[ ,k]
  
  lp_z <- lpNorm( z, p = 2 )
  
  xrng <- c( lpNorm(c(1, rep(0, n-1)), p = 2), lpNorm(rep(1/n, n), p = 2) )
  plot( x = lp_parts, y = apply(theta_mat, 2, max), xlim = xrng, ylim = range(z) )
  points( x = lp_parts, y = apply(theta_mat, 2, min), xlim = xrng, ylim = range(z) )
  points( x = lp_parts, y = lp_z * lp_parts, col = 3 )
  
  
  range( apply( P, 1, lpNorm, 2 ) )
  range( apply( P, 1, lpNorm, 0.5 ) )
  range( apply( P, 1, lpNorm, 100 ) )
  range( apply( P, 1, entropy, base = 2, exponential = FALSE ) )
  
  
  
  
  ########
  
  # Permutations
  
  idx <- which(lp_parts == rev(sort(lp_parts))[2])
  parts_mat[ ,idx]
  P <- matrix( NA, nrow = N, ncol = nrow(parts_mat) )
  for ( j in 1:N ) {
    P[j, ] <- sample( x = parts_mat[ ,idx], replace = FALSE )
  } 
  P_unique <- P[!duplicated(P), ]
  theta <- as.numeric( (P / n) %*% z )
  theta_unique <- as.numeric( (P_unique / n) %*% z )
  
  s <- sum( (theta - mean(theta))^2 ) / length(theta)
  s_u <- sum( (theta_unique - mean(theta_unique))^2 ) / length(theta_unique)
  s
  s_u
  var( theta )
  var( theta_unique )
  

  lp_parts[idx]^2 * var(z) - 1/n * var(z)
  s_u
  var(theta_unique)
  
  mean(theta_unique)
  mean(z)   # same same
  
  
  
  
  
  
  
  ########
  
  
  
  # boxplot( as.data.frame(theta_mat) ) 
  # 
  # ldens <- apply( theta_mat, 2, density )
  # # ldens <- apply( scale(theta_mat, center = rep(mean(z), ncol(theta_mat)), TRUE), 2, density )
  # plot.ldensity( ldens, fillin = FALSE )
  # lines( density( z) )
  
  
  mu_hat <- apply( theta_mat, 2, mean )
  sds_hat <- apply( theta_mat, 2, sd ) 
  sk_hat <- apply( theta_mat, 2, skewness )
  m3_hat <- apply( theta_mat, 2, function(x) { sum( (x - mean(x))^3 ) / length(x) })
  kurt_hat <- apply( theta_mat, 2, kurtosis )
  
  
  plot( mu_hat )
  abline( h = mean(z), col = 2 ) 
  
  plot( sds_hat )
  abline( h = sd(z), col = 2 ) 
  abline( h = sd(z) / sqrt(n) )
  
  plot( sk_hat )
  abline( h = skewness(z), col = 2 ) 
  
  plot( x = sk_hat, y = m3_hat )
  
  plot( kurt_hat )
  abline( h = kurtosis(z), col = 2 ) 
  

  ###
  plot( x = lp_parts^2 , y = sds_hat^2, ylim = c(0, var(z)) )   # !!!!!!
  points( x = lp_parts^2, y = lp_parts^2 * var(z) - min(lp_parts^2 * var(z)), col = 2 )
  ###
  
  plot( x = ent_parts , y = sds_hat ) 
  plot( x = ent_parts, y = m3_hat )
  plot( x = lp_parts100, y = sk_hat )
  plot( x = lp_parts100, y = m3_hat )
  plot( x = lp_parts3^3, y = m3_hat )
  plot( x = lp_parts5^5, y = m3_hat^2 / sds_hat^(3/2) )
  plot( x = lp_parts5^5, y = m3_hat / sds_hat^3 )
  
  
  
  
  
  plot( x = lp_parts, y = 1 / sk_hat^(3/2) )
  
  
  
  
  
  
  
  sd(z)
  sds_hat[length(sds_hat)]; max(sds_hat)
  lp_parts[1]; lp_parts[length(lp_parts)]
  lp_parts[length(lp_parts)] / lp_parts[1]
  lpNorm( rep(1/n, n), p = 2 )
  min(lp_parts)
  
  
  var( apply( gridnodes, 1, function(w) { sum(w * z) } ) )
  var(z)
  sum( prob_tot * (lp_parts^2 * var(z) - min(lp_parts)) )
  sum( prob_tot * sds_hat^2 )
  
  
  sum( prob_tot * sds_hat )
  sd( apply( gridnodes, 1, function(w) { sum(w * z) } ) )
  sqrt( sum( prob_tot * (lp_parts^2 * var(z) - min(lp_parts)) ) )
  
  
  Boot <- boot::boot( data = z,
                      statistic = meanBoot,
                      R = 10^6 )
  sd(Boot$t); var(Boot$t)
  
  
  plot( lp_parts^2 * var(z) - sds_hat^2 )
  max(sds_hat^2)
  max(lp_parts^2)
  var(z)
  
  
  
  
  plot( x = ent_parts, y = sds_hat )
  plot( x = ent_parts, y = sds_hat^2 )
  
  plot( x = lp_parts, y = sk_hat )
  plot( x = lp_parts^2, y = kurt_hat )
  
  
  
  
  
  xrng <- c( entropy(c(1, rep(0, n-1))), entropy(rep(1/n, n)) )
  plot( x = ent_parts, y = apply(theta_mat, 2, max), xlim = xrng, ylim = range(z) )
  points( x = ent_parts, y = apply(theta_mat, 2, min), xlim = xrng, ylim = range(z) )
  points( x = ent_parts, y = lp_z * lp_parts, col = 3 )
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  
  
  
  
  theta <- apply( gridnodes_unique, 1, function(p) { sum(p * z) } )
  ent <- apply( gridnodes_unique, 1, entropy )
  # ent <- apply( gridnodes_unique, 1, entropy, base = 2 )
  # ent <- apply( gridnodes_unique, 1, entropy, exponential = FALSE )
  # ent <- apply( gridnodes_unique, 1, entropy, exponential = FALSE, base = 2 )
  lp <- apply( gridnodes_unique, 1, lpNorm, p = 2 )
  
  
  plot( x = ent, y = theta, ylim = range(z) )
  arrows( x0 = min(ent), y0 = max(z), x1 = max(ent), y1 = mean(z), col = 2 )
  arrows( x0 = min(ent), y0 = min(z), x1 = max(ent), y1 = mean(z), col = 2 )
  
  plot( x = lp, y = theta, ylim = range(z) )
  
  
  
  ###
  xrng <- c( entropy(c(1, rep(0, n-1))), entropy(rep(1/n, n)) )
  plot( x = ent, y = theta, ylim = range(z), xlim = xrng )
  arrows( x0 = 1, y0 = max(z), x1 = n, y1 = mean(z), col = 2 )
  arrows( x0 = 1, y0 = min(z), x1 = n, y1 = mean(z), col = 2 )
  ###
  
  ###
  xrng <- c( lpNorm(c(1, rep(0, n-1)), p = 2), lpNorm(rep(1/n, n), p = 2) )
  plot( x = lp, y = theta, ylim = range(z), xlim = xrng )
  arrows( x0 = max(xrng), y0 = max(z), x1 = min(xrng), y1 = mean(z), col = 2 )
  arrows( x0 = max(xrng), y0 = min(z), x1 = min(xrng), y1 = mean(z), col = 2 )
  ###
  
  
  
  plot( density(z) )
  
  
  
  lTheta <- list()
  for ( i in seq(along = id_match_unique) ) {
    idx <- which(id_match == id_match_unique[i])
    lTheta[[i]] <- apply( matrix(gridnodes_unique[idx, ], nrow = length(idx)), 
                          1, function(p) { sum(p * z) } )
  }
  length(lTheta)
  
  idx <- unlist( lapply( lTheta, FUN = function(x) { length(x) > 1 } ) )
  ldens <- lapply( lTheta[idx], FUN = density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  dim(parts_mat)
  
  
 
  
  
  
  # Cauchy-Schwarz inequality per level-set
  lp_z <- lpNorm(z, p = 2)
  csi <- lp_z * lp_parts[id_match_unique]
  
  plot( x = csi, y = unlist( lapply( lTheta, max) ), xlim = range(z), ylim = range(z) )
  abline( a = 0, b = 1 )
  
  
  
  xrng <- c( lpNorm(c(1, rep(0, n-1)), p = 2), lpNorm(rep(1/n, n), p = 2) )
  plot( x = lp_parts[id_match_unique], y = unlist( lapply( lTheta, min) ),
        xlim = xrng, ylim = range(z) )
  points( x = lp_parts[id_match_unique], y = unlist( lapply( lTheta, max) ),
          xlim = xrng, ylim = range(z) )
  points( x = 1, y = min(z) )
  points( x = 1, y = max(z) )
  abline( h = mean(z) )
  range( unlist( lapply( lTheta, mean ) ) )
  points( x = lp_parts[id_match_unique], y = csi, col = 2 )
  points( x = lp_parts, y = lp_z * lp_parts, col = 3 )
  
  
  
  xrng <- c( entropy(c(1, rep(0, n-1))), entropy(rep(1/n, n)) )
  plot( x = ent_parts[id_match_unique], y = unlist( lapply( lTheta, min) ),
        xlim = xrng, ylim = range(z) )
  points( x = ent_parts[id_match_unique], y = unlist( lapply( lTheta, max) ) )
  points( x = 1, y = min(z) )
  points( x = 1, y = max(z) )
  points( x = ent_parts[id_match_unique], y = csi, col = 2 ) 
  points( x = ent_parts, y = lp_z * lp_parts, col = 3 )
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  
  
  length(id_match_unique)
  sum( prob_tot[ id_match_unique ] )
  plot( cumsum(prob_tot[ id_match_unique ]) )
  
  
  
  
  ent_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, entropy )
  lp_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 2 )
  
  
  lTheta <- list()
  for ( i in seq(along = id_match_unique) ) {
    idx <- which(id_match == id_match_unique[i])
    lTheta[[i]] <- apply( matrix(gridnodes_unique[idx, ], nrow = length(idx)), 
                          1, function(p) { sum(p * z) } )
  }
  length(lTheta)
  
  
  lTheta
  z  
  mean(z)
  parts_mat
  lapply( lTheta, FUN = range )
  unlist(lapply( lTheta, FUN = min )) * prob_tot
  sum( unlist(lapply( lTheta, FUN = min )) * prob_tot[id_match_unique] )
  sum( unlist(lapply( lTheta, FUN = max )) * prob_tot[id_match_unique] )
  range(z)
  
  plot( density( unlist(lTheta) ) )
  abline( v = sum( unlist(lapply( lTheta, FUN = mean )) * prob_tot[id_match_unique] ) )
  abline( v =  sum( unlist(lapply( lTheta, FUN = min )) * prob_tot[id_match_unique] ) )
  abline( v = sum( unlist(lapply( lTheta, FUN = max )) * prob_tot[id_match_unique] ) )
  
  
  tmp <- rbind( parts_mat[ ,id_match_unique], 
                unlist(lapply( lTheta, FUN = min )),
                unlist(lapply( lTheta, FUN = max )), 
                ent_parts[id_match_unique],
                prob_tot[id_match_unique] )
  headleft(tmp)
  tailleft(tmp)
  
  
  
  
  plot( x = ent_parts[id_match_unique], 
        y = unlist( lapply(lTheta, sd) ) )
  
  plot( x = ent_parts[id_match_unique], 
        y = lp_parts[id_match_unique] )
  
  plot( x = log(prob_tot[id_match_unique]), 
        y = ent_parts[id_match_unique] )
  
  
  plot( lp_parts[1:10^3] )
  
  plot( density( lp_parts ) )    
  lines( density( lp), col = 2 )  
  
  plot( density( ent_parts ) )    
  lines( density( ent), col = 2 )  
  
  
  
  
  
  # Cauchy-Schwarz inequality
  
  lp_z <- lpNorm( x = z, p = 2 )
  cbind( unlist( lapply( lTheta, FUN = function(x) { max(x^2) } ) ),
         lp_parts[id_match_unique]^2 * lp_z^2 )
  
  csi <- sum( lp_parts[id_match_unique]^2 * lp_z^2 * prob_tot[id_match_unique] )
  csi
  
  plot( density(z) )
  abline( v = max(z) )
  abline( v = sqrt(csi), col = 3 )
  abline( v = csi, col = 2 )  
  
  plot( density(unlist(lTheta)) )
  lines( density(theta), col = 2 )
  abline( v = sqrt(csi) )  
  abline( v = sum( lpNorm(rep(1/n, n), p = 2) * lp_z ) )
  
  
  sum( lpNorm(rep(1/n, n), p = 2)^2 * lp_z^2 )
  sqrt( sum( lpNorm(rep(1/n, n), p = 2)^2 * lp_z^2 ) )
  sum( lpNorm(rep(1/n, n), p = 2) * lp_z )
  mean(z)  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Check  ||\rho/n||_2^2 * var(Z)
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(GPO)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  lpNorm <- function(x, x_bar = rep(0, length(x)), p = 1) { (sum(abs(x - x_bar)^p))^(1/p) }
  
  
  n <- 3
  N <- 10^4
  set.seed(1111)
  x <- rnorm(n)
  z <- x^2
 
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]

  # Entropy and lp-norm
  ent_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, entropy, exponential = FALSE, base = exp(1) )
  lp_parts <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 2 )
  lp_parts3 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 3 )
  lp_parts5 <- apply( parts_mat[ ,1:min(ncol(parts_mat), 10^5)] / n, 2, lpNorm, p = 5 )
  
  
  # Permutations
  lP <- list()
  for ( k in 1:ncol(parts_mat) ) {
    parts_mat[ ,k]
    P <- matrix( NA, nrow = N, ncol = nrow(parts_mat) )
    for ( j in 1:N ) {
      P[j, ] <- sample( x = parts_mat[ ,k], replace = FALSE )
    } 
    P_unique <- P[!duplicated(P), ]
    lP[[k]] <- P_unique
  }
  
  # Check 
  rbind( draws, lapply(lP, nrow) )
  
  # Compute statistic
  lTheta <- lapply( lP, FUN = function(P) { as.numeric( (P / n) %*% z ) } )
  
  # Variance of the statistcs
  sigma_theta <- unlist( lapply( lTheta, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } ) )
  
  # Skewness of the statistcs
  sk_theta <- unlist( lapply( lTheta, FUN = skewFUN ) )
  
  
  # same same
  lp_parts^2 * var(z) - 1/n * var(z)
  sigma_theta
  # unlist( lapply( lTheta, var ) )
  
  length( unique(lp_parts) )
  length( unique(sigma_theta) )
  
  mean(unlist(lapply(lTheta, mean)))
  mean(z)   # same same
  
  
  plot( x = lp_parts^2, y = lapply(lTheta, var) )
  plot( x = lp_parts^2, y = sigma_theta )
  
  
  m2 <- lp_parts^2 * var(z) - 1/n * var(z)
  m3 <- M3FUN( x = z, M2 = m2 )
  plot( x = m3, y = sk_theta )
  plot( x = m2, y = sigma_theta )
  
  
  cumsum( prob_tot * m3 )
  M3FUN( x = z, M2 = var(z) * ( 2 / (n + 1) - 1 / n ) )
  
  
  
  
  
  idx <- 1
  lp_parts[idx]^2
  lP[[idx]]
  lTheta[[idx]]
  lp_parts[idx]^2 * var(z) - 1/n * var(z)
  sigma_theta[idx]
  
  
  
  apply( lP[[idx]] / n, 2, sum )
  apply( lP[[idx]] / n, 2, mean )
  t( apply( lP[[idx]] / n, 1, function(p) { (p * z)^1 } ) )
  apply( t( apply( lP[[idx]] / n, 1, function(p) { p * z } ) ), 1, sum )
  apply( t( apply( lP[[idx]] / n, 1, function(p) { p * z } ) ), 2, mean )
  z
  
  
  
  ##########
  
  # Step-by-step
  
  idx <- 2
  P <- lP[[idx]] / n
  lM <- lapply( 1:nrow(P), FUN = function(i) { p <- as.numeric(P[i, ]);  matrix(p*z, ncol = 1) %*% matrix(p*z, nrow = 1) })
  lM
  A <- mean( unlist( lapply( lM, sum ) ) )
  B <- (sum(z) / n)^2
  A
  B
  A - B
  theta <- apply( P, 1, function(p) { sum(p * z) } )
  theta
  lTheta[idx]
  v_theta <- sum( (theta - mean(theta))^2 ) / length(theta)
  sum( theta^2 ) / length(theta) - mean(theta)^2 
  v_theta
  sigma_theta[idx]
  
  lapply( 1:nrow(P), FUN = function(p) { sum(p * z) - mean(z) } )
  unlist( lapply( 1:nrow(P), FUN = function(i) { p <- as.numeric(P[i, ]); sum(p * z) - mean(z) } ) )
  mean( unlist( lapply( 1:nrow(P), FUN = function(i) { p <- as.numeric(P[i, ]); sum(p * z) - mean(z) } ) )^2 )
  
  
  
  
  # -------
  
  rho <- parts_mat[ ,idx] / n
  sum(rho^2)
  lpNorm( rho, p = 2 )^2 * var(z)
  sum(rho^2) * 1 / (n - 1) * sum((z - mean(z))^2)
  # sum(rho^2) * 1 / (n - 1) * sum(z^2) - sum(rho^2) * (sum(z) / (1 - n))^2 # wrong!
  sum(rho^2) * (1 / (n - 1) * sum(z^2) - 2 / (n - 1) * sum(z * mean(z)) + mean(z)^2)
  
  D <- sum(rho^2) * (1 / (n - 1) * sum(z^2) - 2 / (n - 1) * sum(z * mean(z)) + mean(z)^2 * n / (n-1))
  E <- (1 / n) * (1 / (n - 1) * sum(z^2) - 2 / (n - 1) * sum(z * mean(z)) +  mean(z)^2 * n / (n-1))
  1 / (n * (n - 1) ) * (sum(z^2) - 2 * sum(z * mean(z)) + n * mean(z)^2)
  D
  E
  
  A - B
  D - E
  
  A; D
  B; E
  A - D
  B - E
  
  
  p <- parts_mat[ ,idx] / n
  var(z) * (sum(p^2) - 1/n)
  1 / nrow(Q) * sum( (theta - mean(theta))^2 )  # same same
  
  sum( (z - mean(z))^2 ) 
  sum( (theta - mean(theta))^2 )
  
  scl <- (sum(p^2) - 1/n) / (n - 1)
  sum( (z - mean(z))^2 ) * scl
  sum( (theta - mean(theta))^2 ) / nrow(P)
  
  a <- unlist( lapply( 1:nrow(P), FUN = function(i) { p <- as.numeric(P[i, ]); sum(p * z) } ) )
  mean( (a - mean(z))^2 )
  
  
  a
  z
  P
  p
  
  tmp <- c( sum(z^2), 2 * mean(z) * sum(z), n * mean(z)^2 )
  tmp
  tmp * scl
  tmp[1] - tmp[2] + tmp[3]
  (tmp[1] - tmp[2] + tmp[3]) * scl
  
  
  
  
  
  
  ###
  
  # Toy example
  
  n <- 3
  z <- c(98, 100, 1111)
  p <- c(2, 1, 0)
  Q <- rbind( c(1, 2, 0), 
              c(0, 1, 2), 
              c(2, 0, 1), 
              c(1, 0, 2), 
              c(2, 1, 0), 
              c(0, 2, 1) )
  p <- p / n
  Q <- Q / n
  theta <- apply( Q, 1, function(q) { sum(q * z) } )
  
  
  var(z) * (sum(p^2) - 1/n)
  1 / nrow(Q) * sum( (theta - mean(theta))^2 )  # same same
 
  sum( (z - mean(z))^2 ) * (sum(p^2) - 1/n) / (n - 1)
  sum( (theta - mean(theta))^2 ) / nrow(Q)
  
  
  
  parts_mat <- partitions::parts(n = n)
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  
  
  
  # Skewness
  m2 <- var(z) * (sum(p^2) - 1/n)
  m3 <- M3FUN( x = z, M2 = m2 )
  m3
  skewFUN(theta)
    
  
  
  
  
  
