  
  
  ############################################################################
  ### EXACT MOMENTS - INFORMATIVE PRIOR - ASYMMETRY
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     15.05.2022
  # First version:    15.05.2022
  # --------------------------------------------------------------------------
  
  
  # Sections:
  #
  # Bayesian bootstrap
  # Face-Edge-Vertex (FEV) bias
  
  
  
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
  
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  
  
  
  
  betaFUN <- function(a, a0, m) { betaMoments( a = a, a0 = a0, m = m ) }
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # Loop over n and scale factor
  # --------------------------------------------------------------------------
  
  
  # n_vec <- round( seq(from = 2, to = 10^3, length.out = 100) )
  # n_vec <- 2:10^3
  n_vec <- c(2:10, 20, 50, 100, 200, 500, 1000)
  # lambda_vec <- c( seq(from = 0.001, to = 1, length.out = 100/2), 
  #                  (1 + log(seq(from = 1.001, to = 10^3, length.out = 100/2)))^2 )
  lambda_vec <- exp( seq(from = -10, to = 10, length.out = 100) )
  
  EE2N2 <- matrix( NA, length(n_vec), length(lambda_vec) )
  EE3N3 <- EE2N2
  EE4N4 <- EE2N2
  EE2N2Wbar <- EE2N2 
  
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    alpha_base <- rep(1, n)
    for ( j in seq(along = lambda_vec) ) {
      lambda <- lambda_vec[j]
      alpha <- alpha_base * lambda
      # Exact expected p-norm (to the power of p)
      EE2N2[i, j] <- sum( unlist( lapply( alpha, FUN = betaFUN, a0 = sum(alpha), m = 2 ) ) )
      EE3N3[i, j] <- sum( unlist( lapply( alpha, FUN = betaFUN, a0 = sum(alpha), m = 3 ) ) )
      EE4N4[i, j] <- sum( unlist( lapply( alpha, FUN = betaFUN, a0 = sum(alpha), m = 4 ) ) )
      
      # EE2N2Wbar[i, j] <- sum(rep(1/n, n)^2)
      EE2N2Wbar[i, j] <- 1/n
    }  
  }
  
  headleft( EE2N2 )
  headleft( EE2N2Wbar )
  headleft( EE2N2 - EE2N2Wbar )
  
  # require(rgl)
  # plot3d( x = n_vec, y = log(lambda_vec), z = EE2N2 )
  # plot3d( x = 1:length(n_vec), y = 1:length(lambda_vec), z = EE2N2 )
  
  
  
  
  
  
  colors <- fBasics::divPalette(n = nrow(EE2N2), "RdYlGn")
  plot( x = log(lambda_vec), y = EE2N2[1, ], ylab = "E(||w||_2^2)", xlab = "log(lambda)", type = "l", ylim = c(0, 1) )
  abline( v = 0, col = "grey" )
  abline( h = 0.5, col = "grey" )
  for ( i in 1:nrow(EE2N2) ) {
    points( x = log(lambda_vec), y = EE2N2[i, ], col = colors[i], type = "l" )
  }
  legend("topright", lwd= 2, legend = paste0("n = ", n_vec), col = colors, text.col = colors, bty = "n", cex = 0.7)
 
  
  
  plot( x = log(lambda_vec), y = EE3N3[1, ], ylab = "E(||w||_3^3)", xlab = "log(lambda)", 
        type = "l", ylim = c(0, 1) )
  abline( v = 0, col = "grey" )
  abline( h = 0.5, col = "grey" )
  for ( i in 1:nrow(EE3N3) ) {
    points( x = log(lambda_vec), y = EE3N3[i, ], col = colors[i], type = "l" )
  }
  
  
  
  plot( x = log(lambda_vec), y = EE4N4[1, ], ylab = "E(||w||_4^4)", xlab = "log(lambda)", 
        type = "l", ylim = c(0, 1) )
  abline( v = 0, col = "grey" )
  abline( h = 0.5, col = "grey" )
  for ( i in 1:nrow(EE4N4) ) {
    points( x = log(lambda_vec), y = EE4N4[i, ], col = colors[i], type = "l" )
  }

    
  
  
  
  colors <- fBasics::divPalette(n = nrow(EE2N2), "RdYlGn")
  plot( x = log(lambda_vec), y = EE2N2[1, ] - EE2N2Wbar[1, ], type = "l", ylim = c(0, 1) )
  abline( v = 0, col = "grey" )
  abline( h = 0.5, col = "grey" )
  for ( i in 1:nrow(EE2N2) ) {
    points( x = log(lambda_vec), y = EE2N2[i, ] - EE2N2Wbar[i, ], col = colors[i], type = "l" )
  }
  
  xxx <- log(lambda_vec)
  tmp <- 1 / (1 + exp(xxx)) / 2
  points( x = xxx, y = tmp )  
  
  
  
  
  
  # Asymmetric \alpha, sum(a) = n
  
  n_vec <- c(2:10, 20, 50, 100, 200, 500, 1000)
  fev_vec <- 1:20
  EE2N2 <- matrix( NA, length(n_vec), length(fev_vec) )
  EE2N2Wbar <- EE2N2
  
  set.seed(1234)
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    # alpha_base <- as.numeric( rdirichlet(n = 1, alpha = runif(n)) )
    alpha_base <- rep(1/n, n)
    eps <- runif(n) / 20
    for ( j in seq(along = fev_vec) )  {
      alpha <- fevBias( x = alpha_base + eps, q = fev_vec[j] ) * n
      # Exact expected p-norm (to the power of p)
      betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
      beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
      EE2N2[i, j] <- sum(beta_moments)
      EE2N2Wbar[i, j] <- sum( (alpha/n)^2 )
    }
  }
  
  
  colors <- fBasics::divPalette(n = nrow(EE2N2), "RdYlGn")
  plot( x = fev_vec, y = EE2N2[1, ] - EE2N2Wbar[1, ], type = "l", ylim = c(0, 0.5) )
  for ( i in 1:nrow(EE2N2) ) {
    points( x = fev_vec, y = EE2N2[i, ] - EE2N2Wbar[i, ], col = colors[i], type = "l" )
  }
  
  
  
  
  
  
  # Asymmetric \alpha, sum(a) <> n
  
  n_vec <- c(2:10, 20, 50, 100, 200, 500, 1000)
  lambda_vec <- exp( seq(from = -10, to = 10, length.out = 100) )
  EE2N2 <- matrix( NA, length(n_vec), length(lambda_vec) )
  EE2N2Wbar <- EE2N2
  relent <- rep(NA, length(n_vec))
  fev_vec <- 1:10
  
  set.seed(1234)
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    alpha_base <- as.numeric( rdirichlet(n = 1, alpha = runif(n)) )
    relent[i] <- entropy( p = alpha_base ) / n
    for ( k in seq(along = fev_vec) )  {
      alpha_base_k <- fevBias( x = alpha_base, q = fev_vec[k] )
    }
    for ( j in seq(along = lambda_vec) ) {
      lambda <- lambda_vec[j]
      alpha <- alpha_base_k * lambda
      # Exact expected p-norm (to the power of p)
      betaFUN <- function(a) { betaMoments( a = a, b = sum(alpha) - a, m = 2 ) }
      beta_moments <- unlist( lapply( alpha, FUN = betaFUN ) )
      EE2N2[i, j] <- sum(beta_moments)
      EE2N2Wbar[i, j] <- sum((alpha_base/n)^2)
    }  
  }
  
  headleft( EE2N2 )
  headleft( EE2N2Wbar )
  headleft( EE2N2 - EE2N2Wbar )
  
  
  
  colors <- fBasics::divPalette(n = nrow(EE2N2), "RdYlGn")
  plot( x = log(lambda_vec), y = EE2N2[1, ] - EE2N2Wbar[1, ], type = "l", ylim = c(0, 1) )
  abline( v = 0, col = "grey" )
  abline( h = 0.5, col = "grey" )
  for ( i in 1:nrow(EE2N2) ) {
    points( x = log(lambda_vec), y = EE2N2[i, ] - EE2N2Wbar[i, ], col = colors[i], type = "l" )
  }
  
  
  # plot( x = relent, y = EE2N2[ ,1] - EE2N2Wbar[ ,1], type = "l", ylim = c(0, 1) )
  # for ( i in 1:ncol(EE2N2) ) {
  #   points( x = relent, y = EE2N2[ ,i] - EE2N2Wbar[ ,i], col = colors[i], type = "l" )
  # }
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Face-Edge-Vertex (FEV) bias
  # --------------------------------------------------------------------------
  
  
  n_vec <- c(2:10, 20, 50, 100, 200, 500, 1000)
  # q_vec <- seq(from = -10, to = 10, length.out = 100)
  q_vec <- seq(from = 0, to = 100, by = 1)
  # n_vec <- round(seq(from = 2, to = 100, length.out = 10))
  # q_vec <- c(1, seq(from = 10, to = 100, by = 10) )
  n_sim <- 10^3
  EE2N2 <- matrix( NA, n_sim, length(q_vec) )
  EE2N2_ana <- EE2N2
  ee2n2 <- numeric(length(n_vec) )
  lEE2N2 <- list()
  lEE2N2_ana <- list()
  
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i]
    alpha <- rep(1, n)
    ee2n2[i] <- sum( unlist( lapply( alpha, FUN = betaFUN, m = 2 ) ) )
    P <- rdirichlet( n = n_sim, alpha = alpha )
    for ( j in seq(along = q_vec) ) {
      P_fev <- fevBias( x = P, q = q_vec[j] )
      EE2N2[ ,j] <- apply( P_fev, 1, function(p) { sum(p^2) } )
      num <- unlist( lapply( alpha, FUN = betaFUN, m = 2 * q_vec[j] ) )
      denom <- sum( unlist( lapply( alpha, FUN = betaFUN, m = q_vec[j] ) ) )^2
      EE2N2_ana[ ,j] <- sum(num / denom)
    }
    lEE2N2[[i]] <- EE2N2
    lEE2N2_ana[[i]] <- EE2N2_ana
  }
  
  tmp <- lapply( lEE2N2, function(A) { apply(A, 2, mean) } )
  ee2n2_mat <- do.call( cbind, tmp )
  dim(ee2n2_mat)
  
  
  
  
  n <- 10
  alpha <- rep(1, n)
  P <- rdirichlet( n = 10^5, alpha = alpha )
  mean( apply( P, 1, function(p) { sum(p^2) } ) )
  sum( unlist( lapply( alpha, FUN = betaFUN, m = 2 ) ) )
  
  
  P_fev <- fevBias( x = P, q = 4 )
  num <- unlist( lapply( alpha, FUN = betaFUN, m = 2 * 4 ) )
  denom <- sum( unlist( lapply( alpha, FUN = betaFUN, m = 4 ) ) )^2
  sum(num / denom)
  mean( apply( P_fev, 1, function(p) { sum(p^2) } ) )
  
  
  
  
  
  
  
  colors <- fBasics::divPalette(n = ncol(ee2n2_mat) )
  plot( log(q_vec), ee2n2_mat[ ,1], ylim = c(0, 1), type = "l" )
  for ( j in 1:ncol(ee2n2_mat) ) {
    lines( log(q_vec), ee2n2_mat[ ,j], col = colors[j] )
  }
  points( x = rep(0, length(ee2n2)), y = ee2n2, col = colors, pch = 19 )
  
  
  
  
  plot3d( x = q_vec, y = n_vec, z = ee2n2_mat, col = colors )
  
  surface3d( x = q_vec * 100, y = n_vec, z = ee2n2_mat * 1000, col = colors )
  
  
  
  tmp <- lapply( lEE2N2_ana, function(A) { apply(A, 2, mean) } )
  ee2n2_ana_mat <- do.call( cbind, tmp )
  
  headleft(ee2n2_mat)
  headleft(ee2n2_ana_mat)
  
  
  
  
  colors <- fBasics::divPalette(n = ncol(ee2n2_mat) )
  plot( q_vec, ee2n2_ana_mat[ ,1], ylim = c(0, max(na.omit(ee2n2_ana_mat))), type = "l" )
  plot( q_vec, ee2n2_ana_mat[ ,1], ylim = c(0, 10), type = "l" )
  for ( j in 1:ncol(ee2n2_ana_mat) ) {
    lines( q_vec, ee2n2_ana_mat[ ,j], col = colors[j] )
  }
  points( x = rep(1, length(ee2n2)), y = ee2n2, col = colors, pch = 19 )
  
  
  
  
  
  # What is the limes n --> \Infty E(||\omega||_r^r) ?
  
  n_vec <- 1:100
  ans <- matrix(NA, nrow = length(n_vec), ncol = 5)
  for ( i in seq(along = n_vec) ) {
    n <- n_vec[i] 
    ans[i, 1] <- n * betaMoments(a = 1, b = n-1, m = 2)
    ans[i, 2] <- n * betaMoments(a = 1, b = n-1, m = 3)
    ans[i, 3] <- n * betaMoments(a = 1, b = n-1, m = 4)
    ans[i, 4] <- n * betaMoments(a = 1, b = n-1, m = 5)
    ans[i, 5] <- n * betaMoments(a = 1, b = n-1, m = 6)
  }
  
  # colors <- fBasics::divPalette(n = ncol(ans))
  colors <- 1:5
  plot( x = log(n_vec), y = ans[ ,1], type = "o" )
  for ( j in 1:ncol(ans) ) {
    lines( x = log(n_vec), y = ans[ ,j], type = "o", col = colors[j] )
  }
  n <- 10^6
  abline( h = n * betaMoments(a = 1, b = n-1, m = 2) )

 
  
  
    
  
  
  
  
  
  