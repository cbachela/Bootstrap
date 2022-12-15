  
  
  ############################################################################
  ### REGRESSION BOOTSTRAP - PAIRS BOOTSTRAP AND WEIGHTED LEAST SQUARES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.07.2022
  # First version:    01.07.2022
  # --------------------------------------------------------------------------
  

  
  
  
  # --------------------------------------------------------------------------
  # Require
  # --------------------------------------------------------------------------
  
  require(partitions)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
  source( paste0(wd, "Regression/helper_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  X_bm <- getMSCIData( universe = "bm", frqncy = "m", ccy = "usd" )
  X <- getMSCIData( universe = "dm", frqncy = "m", ccy = "usd" )
  
  # Fama-French factors
  path_ff <- "R:/Asset_Management/R/myRData/Factor/"
  env_ff <- readRDS( file = paste0(path_ff, "data.Rds") )
  FF3 <- env_ff$lFactor$`3F`$DM_monthly
  FF5 <- env_ff$lFactor$`5F`$DM_monthly
  
  
  TD <- trainingData( Y_train = head(X[ ,1], 10),
                      X_train = X_bm )
  
 
  
  
  
  
  # Synthetic data - normal distribution
  n <- 10^1
  set.seed(1111)
  y <- scale( rnorm( n, 0, 1 ), TRUE, TRUE )
  r <- scale( rnorm( n, 0, 1 ), TRUE, TRUE )
  a <- 0.01
  b <- 0.2
  x <- (y - a - r) / b
  
  TD <- trainingData( Y_train = asTimeSeries(matrix(y, ncol = 1)),
                      X_train = asTimeSeries(matrix(x, ncol = 1)) )
  
  # Checks
  mean(y)
  mean(r)
  mean(x)
  cbind(y, a + b*x + r)
  
  
  
  # --------------------------------------------------------------------------
  # OLS Regression
  # --------------------------------------------------------------------------
  
  Y_train <- TD$Y_train[ ,1]
  X_train <- TD$X_train[ ,1]
  reg <- regression( Y_train = Y_train,
                     X_train = X_train,
                     type = "ols" )
  coeff <- reg$coeffmat
  a_hat <- coeff[1, 1]
  b_hat <- coeff[2, 1]
  resid <- residuals(reg$reg)
  summary(reg$reg)  
  
  
  # OLS
  betaHatUnivariate( x = TD$X_train, y = TD$Y_train, 
                     x_bar = mean(TD$X_train), y_bar = mean(TD$Y_train) )
  betaHat( X = TD$X_train, y = TD$Y_train )
  
  # WLS
  n <- nrow(TD$X_train)
  wghts <- rep(1, n)
  wghts <- (1:n) / sum(1:n)
  betaHat( X = TD$X_train, y = TD$Y_train, wghts = wghts, type = "wls" )
  
  
  
  
  # --------------------------------------------------------------------------
  # Bootstrap (classical) of regression model by resampling Z_t = (Y, X)_t
  # --------------------------------------------------------------------------
  
  n_sim <- 10^3
  n <- nrow(TD$Y_train)
  gridnodes <- matrix( 0, nrow = n_sim, ncol = n )
  ab_boot <- matrix( 0, nrow = n_sim, ncol = 2, dimnames = list(NULL, c("a", "b")) )
  ab_boot2 <- ab_boot
  ab_ols <- ab_boot[ ,rep(1, 4)]
  ab_wls <- ab_boot[ ,rep(1, 4)]
  colnames(ab_ols) <- colnames(ab_wls) <- c("a", "b", "b_num", "b_denom")
  
  for ( i in 1:n_sim ) {
    
    # Boot
    idx <- sort( sample( x = 1:n, size = n, replace = TRUE ) )
    #// Create synthetic (future) time stamps bcs dates need to be unique
    Y_train_boot <- matrix( TD$DF[idx, 1], ncol = 1 )
    X_train_boot <- matrix( TD$DF[idx, -1], ncol = ncol(TD$DF)-1 )
    DF_boot <- asTimeSeries( cbind(Y_train_boot, X_train_boot) )
    reg_tmp <- regression( Y_train = DF_boot[ ,1],
                           X_train = DF_boot[ ,-1],
                           type = "ols" )
    ab_boot[i, ] <- reg_tmp$coeffmat[1:2, 1]
  
    # OLS
    # debugonce( betaHat )
    ab_ols[i, c("b", "b_num", "b_denom")] <- betaHat( X = reg_tmp$EV,
                                                      y = reg_tmp$RV )
    ab_ols[i, "a"] <- alphaHat( y_bar = mean(reg_tmp$RV), 
                                b_hat = ab_ols[i, "b"], 
                                x_bar = mean(reg_tmp$EV) )
    
    # WLS
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    wghts <- gridnodes[i, ]
    # debugonce( betaHat )
    ab_wls[i, c("b", "b_num", "b_denom")] <- betaHat( X = TD$X_train,
                                                      y = TD$Y_train,
                                                      wghts = wghts,
                                                      type = "wls" )
    y_bar <- sum( wghts * TD$Y_train )
    x_bar <- sum( wghts * TD$X_train )
    ab_wls[i, "a"] <- alphaHat( y_bar = y_bar, 
                                b_hat = ab_wls[i, "b"], 
                                x_bar = x_bar )
    
    # This should give different results
    DF_boot2 <- as.timeSeries( apply( TD$DF, 2, function(x) { x * wghts } ) )
    reg_tmp2 <- regression( Y_train = DF_boot2[ ,1],
                            X_train = DF_boot2[ ,-1],
                            type = "ols" )
    ab_boot2[i, ] <- reg_tmp2$coeffmat[1:2, 1]
    
  }
  
  
  a_mat <- cbind( boot = ab_boot[ ,"a"], ols = ab_ols[ ,"a"], wls = ab_wls[ ,"a"], boot2 = ab_boot2[ ,"a"] )
  b_mat <- cbind( boot = ab_boot[ ,"b"], ols = ab_ols[ ,"b"], wls = ab_wls[ ,"b"], boot2 = ab_boot2[ ,"b"] )
  
  
  coeff[1, ]
  apply( a_mat, 2, mean ); apply( a_mat, 2, sd )
  barplot( abs( coeff[1, 1] - apply( a_mat, 2, mean ) ) )
  barplot( abs( coeff[1, 2] - apply( a_mat, 2, sd ) ) )
  
  coeff[2, ]
  apply( b_mat, 2, mean ); apply( b_mat, 2, sd )
  barplot( abs( coeff[2, 1] - apply( b_mat, 2, mean ) ) )
  barplot( abs( coeff[2, 2] - apply( b_mat, 2, sd ) ) )
  
  
  
  ldens_a <- apply( a_mat, 2, density )
  ldens_b <- apply( b_mat, 2, density )
  
  slolz:::plot.ldensity( ldens_a, fillin = FALSE, col = 1:length(ldens_a) )
  abline( v = coeff[1, 1] )

  slolz:::plot.ldensity( ldens_b, fillin = FALSE, col = 1:length(ldens_b) )
  abline( v = coeff[2, 1] )
  abline( v = betaHat( X = TD$X_train, y = TD$Y_train ), col = 4 )
  
  
  
  
  ab_ols
  
  ldens <- apply( ab_ols, 2, densFUN )
  slolz:::plot.ldensity( ldens[-1], fillin = FALSE, col = 1:3 )    

  
  
  
  
  x <- as.numeric(TD$X_train)
  y <- as.numeric(TD$Y_train)
  n <- length(x)
  apply( cbind(x*y, x^2), 2, var) * ( 2/n - 1/n^2 - 1/n )
  apply( ab_ols[ ,-1], 2, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } )
  apply( ab_wls[ ,-1], 2, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } )
  
  
  M2Exact( z = x*y, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = x^2, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = (x*y) / x^2, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = NULL, m2 = var(x*y) / var(x^2), w_bar = rep(1/n, n), exp2norm2 = 2/n - 1/n^2 )
  
  
  M2Exact( z = x*y, exp2norm2 = 2/n - 1/n^2 )
  var( ab_ols[ ,"b_num"] )
  
  
  # Why is  
  # M2Exact( z = x*y, exp2norm2 = 2/n - 1/n^2 )
  # not the same as
  # var( ab_ols[ ,"b_num"] )

  
  
  
  
  n_sim <- 10^3
  n <- nrow(TD$Y_train)
  b_ols <- matrix( 0, nrow = n_sim, ncol = 3, dimnames = list(NULL, c("b", "b_num", "b_denom") ) )
  b_olsU <- b_ols
  
  for ( i in 1:n_sim ) {
    
    idx <- sort( sample( x = 1:n, size = n, replace = TRUE ) )
    Y_train_boot <- matrix( TD$DF[idx, 1], ncol = 1 )
    X_train_boot <- matrix( TD$DF[idx, -1], ncol = ncol(TD$DF)-1 )
    b_ols[i, c("b", "b_num", "b_denom")] <- betaHat( X = X_train_boot,
                                                     y = Y_train_boot )
    b_olsU[i, c("b", "b_num", "b_denom")] <- betaHatUnivariate( x = X_train_boot,
                                                                y = Y_train_boot )
  }
  
  apply( b_ols, 2, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } )
  apply( b_ols, 2, FUN = var )
  apply( b_olsU, 2, FUN = var )
  
 
  x <- as.numeric(TD$X_train)
  y <- as.numeric(TD$Y_train)
  n <- length(x)
  M2Exact( z = NULL, m2 = var(x*y) / var(x^2), w_bar = rep(1/n, n), exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = x*y, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = x^2, exp2norm2 = 2/n - 1/n^2 )
  M2Exact( z = (x - mean(x))^2, exp2norm2 = 2/n - 1/n^2 )
  M2Analytic( z = x^2, method = "classic" )
  
  z <- (x - mean(x))^2
  M2Exact( z = z, exp2norm2 = 2/n - 1/n^2 )
  
  M3Exact( )
  
  
  
  # Bayesian 
  
  n_sim <- 10^3
  n <- nrow(TD$DF)
  samples <- rdirichlet( n = n_sim, alpha = rep(1, n) )
  b_ols <- matrix( 0, nrow = n_sim, ncol = 3, dimnames = list(NULL, c("b", "b_num", "b_denom") ) )
  b_olsU <- b_ols
  theta <- b_ols[ ,1]
  
  for ( i in 1:n_sim ) {
    
    b_ols[i, c("b", "b_num", "b_denom")] <- betaHat( X = TD$X_train,
                                                     y = TD$Y_train,
                                                     wghts = as.numeric(samples[i, ]),
                                                     type = "wls" )
    b_olsU[i, c("b", "b_num", "b_denom")] <- betaHatUnivariate( x = TD$X_train,
                                                                y = TD$Y_train,
                                                                wghts = as.numeric(samples[i, ]),
                                                                type = "wls" )
  }
  x <- as.numeric(TD$X_train)
  y <- as.numeric(TD$Y_train)
  z <- (x - mean(x)) * (y - mean(y)) / (x - mean(x))^2
  theta <- apply( samples, 1, function(w) { sum(w * z) } )
  var(theta)
  
  apply( b_ols, 2, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } )
  
  
  M2Exact( z = NULL, m2 =  var( (x - mean(x)) * (y - mean(y)) ) / var( (x - mean(x))^2 ), w_bar = rep(1/n, n), exp2norm2 = 2 / (n + 1) )
  M2Exact( z = (x - mean(x)) * (y - mean(y)), exp2norm2 = 2 / (n +1) )
  M2Exact( z = x^2, exp2norm2 = 2 / (n +1) )
  M2Exact( z = (x - mean(x))^2, exp2norm2 = 2 / (n +1) )
  
  
  mean( apply( samples, 1, function(x) { sum(x^2) } ) )
  2 / (n + 1)
  
  
  M2Exact( z = NULL, m2 =  var( (x - mean(x)) * (y - mean(y)) / (x - mean(x))^2 ), w_bar = rep(1/n, n), exp2norm2 = 2 / (n + 1) )
  M2Exact( z = NULL, m2 =  var( (y - mean(y)) / (x - mean(x)) ), w_bar = rep(1/n, n), exp2norm2 = 2 / (n + 1) )
  
  
  
  var( (x - mean(x)) * (y - mean(y)) ) / var( (x - mean(x))^2 )
  apply( b_ols, 2, FUN = function(x) { sum( (x - mean(x))^2 ) / length(x) } )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Varsi
  # --------------------------------------------------------------------------
  
  x <- as.numeric(TD$X_train)
  y <- as.numeric(TD$Y_train)
  
  q_vec <- seq( from = 0, to = 1, length.out = 10^4 )
  z_vec_x2 <- quantile( (x - mean(x))^2, q_vec )
  z_vec_xy <- quantile( (x - mean(x)) * (y - mean(y)), q_vec )
  
  FUN <- function(z0, z) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) }
  p_vec_x2 <- unlist( lapply( z_vec_x2, FUN, z = (x - mean(x))^2 ) )
  p_vec_xy <- unlist( lapply( z_vec_xy, FUN, z = (x - mean(x)) * (y - mean(y)) ) )
  
  plot( x = z_vec_x2, y = p_vec_x2, type = "o" )
  lines( x = z_vec_xy, y = p_vec_xy, type = "o", col = 2 )
  plot( x = q_vec, y = p_vec_x2, type  = "o" )
  lines( x = q_vec, y = p_vec_xy, type = "o", col = 2 )  
  
  p_vec_ratio <- p_vec_xy / p_vec_x2
  plot( x = q_vec, y = p_vec_ratio, type  = "o" )
  
  
  tmp <- cbind( p_vec_x2, p_vec_xy, p_vec_ratio )
  tmp[5000:5001, ]
  head(tmp)
  
  
  
  
  
  
  # Derivative
  
  idx <- 3:length(p_vec_x2)
  d_vec_x2 <- (p_vec_x2[idx] - p_vec_x2[idx-2]) / abs(z_vec_x2[idx] - z_vec_x2[idx-2])
  d_vec_x2_std <- d_vec_x2 / sum(d_vec_x2)
  
  idx <- 3:length(p_vec_xy)
  d_vec_xy <- (p_vec_xy[idx] - p_vec_xy[idx-2]) / abs(z_vec_xy[idx] - z_vec_xy[idx-2])
  d_vec_xy_std <- d_vec_xy / sum(d_vec_xy)
  
  
  plot( x = z_vec_x2[-c(1, length(z_vec_x2))], y = d_vec_x2_std, type = "o" )
  points( x = z_vec_xy[-c(1, length(z_vec_xy))], y = d_vec_xy_std, type = "o", col = 2 )
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Dot product
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  n <- 10
  M <- 6
  n_boot <- 50
  B <- 10^3 + 1
  
  set.seed(1111)
  df <- 4
  x <- rt( n = n, df = df )
  z <- x^2
  y <- x * rt( n = n, df = 5 )
  
  
  alpha <- rep(1, n)
  w_bar <- alpha / n
  samples <- rdirichlet( n = B, alpha = alpha )
  
  
  wghts <- as.numeric(samples[1, ])
  sum( wghts * z )
  dot(wghts, z)
  sqrt( sum( wghts^2 ) ) * sqrt( sum( z^2 ) ) * cos( angle( v1 = wghts, v2 = z ) )

  theta_z <- apply( samples, 1, angle, v2 = z )  
  theta_y <- apply( samples, 1, angle, v2 = y )
  plot(theta_z)  
  mean(theta_z)
  mean(cos(theta_z))
  cos(mean(theta_z))
  sum( w_bar * z ) / ( sqrt( 2 / (n + 1) ) * sqrt( sum( z^2 ) ) )
  
  mean(cos(theta_y))
  sum( w_bar * y ) / ( sqrt( 2 / (n + 1) ) * sqrt( sum( y^2 ) ) )
  
  
  
  mean( apply( samples, 1, function(x) { sum(x^2) } ) )
  2 / (n + 1)  
  mean( apply( samples, 1, function(x) { sqrt(sum(x^2)) } ) )
  sqrt( 2 / (n + 1) )  
  
  
  
  tmp <- apply( samples, 1, function(w) { sum(w * y) / sum(w * z) } )
  var(tmp)  
  M2Exact( z = NULL, m2 = var( y / z), exp2norm2 = 2 / (n + 1) )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Visualization on simplex
  # --------------------------------------------------------------------------
  
  require(rgl)
  
  n <- 3
  n_sim <- 10^5
  set.seed(1234)
  df <- 5
  # x <- rt( n = n, df = df )
  # y <- rt( n = n, df = df-1 )
  x <- 1:n
  y <- runif(n = n, min = -100, max = 100)
  z1 <- x^2
  z2 <- x * y
  z2z1 <- z2 / z1
  
  
  alpha <- rep(1, n)
  w_bar <- alpha / n
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  theta_z1 <- apply( samples, 1, function(w) { w %*% z1 } )
  theta_z2 <- apply( samples, 1, function(w) { w %*% z2 } )
  theta_ratio <- apply( samples, 1, function(w) { (w %*% z2) / (w %*% z1) } )
  theta_z2z1 <- apply( samples, 1, function(w) { w %*% z2z1 } )
  
  
  
  plot( density(theta_z1) )  
  idx_z1 <- which( abs(theta_z1 - median(theta_z1)) < 1e-03 )
  length(idx_z1)
  
  plot( density(theta_z2) )  
  idx_z2 <- which( abs(theta_z2 - median(theta_z2)) < 1e-03 )
  length(idx_z2)
  
  plot( density(theta_ratio) )  
  idx_ratio <- which( abs(theta_ratio - median(theta_ratio)) < 1e-03 )
  length(idx_ratio)
  
  plot( density(theta_z2z1) )  
  idx_z2z1 <- which( abs(theta_z2z1 - median(theta_z2z1)) < 1e-03 )
  length(idx_z2z1)

  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = "grey", size = 0.5 )
  points3d( samples[idx_z1, 1], y = samples[idx_z1, 2], z = samples[idx_z1, 3], col = "black", size = 3 )
  points3d( samples[idx_z2, 1], y = samples[idx_z2, 2], z = samples[idx_z2, 3], col = "red", size = 3 )
  points3d( samples[idx_ratio, 1], y = samples[idx_ratio, 2], z = samples[idx_ratio, 3], col = "blue", size = 3 )
  points3d( samples[idx_z2z1, 1], y = samples[idx_z2z1, 2], z = samples[idx_z2z1, 3], col = "magenta", size = 3 )
  
  
  
  # Colorcoding
  nq <- 10
  color_vec_z1 <- fBasics::divPalette(n = nq, "RdYlGn" )
  color_vec_z2 <- fBasics::divPalette(n = nq, "RdYlBu" )
  color_vec_ratio <- fBasics::divPalette(n = nq, "Spectral")
  color_vec_z2z1 <- fBasics::divPalette(n = nq, "BrBG")
  
  qvec_z1 <- quantile( theta_z1, seq(from = 0, to = 1, length.out = nq) )
  colors_z1 <- theta_z1 * NA
  for ( i in length(qvec_z1):2 ) {
    colors_z1[ which( theta_z1 <= qvec_z1[i] ) ] <- color_vec_z1[i]
  }
  
  qvec_z2 <- quantile( theta_z2, seq(from = 0, to = 1, length.out = nq) )
  colors_z2 <- theta_z2 * NA
  for ( i in length(qvec_z2):2 ) {
    colors_z2[ which( theta_z2 <= qvec_z2[i] ) ] <- color_vec_z2[i]
  }
  
  qvec_ratio <- quantile( theta_ratio, seq(from = 0, to = 1, length.out = nq) )
  colors_ratio <- theta_ratio * NA
  for ( i in length(qvec_ratio):2 ) {
    colors_ratio[ which( theta_ratio <= qvec_ratio[i] ) ] <- color_vec_ratio[i]
  }

  qvec_z2z1 <- quantile( theta_z2z1, seq(from = 0, to = 1, length.out = nq) )
  colors_z2z1 <- theta_z2z1 * NA
  for ( i in length(qvec_z2z1):2 ) {
    colors_z2z1[ which( theta_z2z1 <= qvec_z2z1[i] ) ] <- color_vec_z2z1[i]
  }
  
  
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = colors_z1, size = 0.5 )
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = colors_z2, size = 0.5 )
  plot3d( x = samples[ ,1], y = samples[ ,2], z = samples[ ,3], col = colors_ratio, size = 0.5 )
  
  
  
  plot( x = samples[ ,1], y = samples[ ,2], col = colors_z1, pch = 19, cex = 0.01 )
  plot( x = samples[ ,1], y = samples[ ,2], col = colors_z2, pch = 19, cex = 0.01 )
  plot( x = samples[ ,1], y = samples[ ,2], col = colors_ratio, pch = 19, cex = 0.01 )
  plot( x = samples[ ,1], y = samples[ ,2], col = colors_z2z1, pch = 19, cex = 0.01 )
  
  
  
  
  
  
  require(RP)
  wd <- "H:/R/RP/R/"
  source( paste0(wd, "class_Simplex.R") )
  source( paste0(wd, "class_Polytope.R") )
  source( paste0(wd, "class_Ellipsoid.R") )
  source( paste0(wd, "geometric_walks.R") )
  
  
  
  
  P <- Polytope$new()
  P$A <- rbind( rep(1, n), diag(rep(1, n)), diag(rep(-1, n)) ) 
  P$b <- c(1, rep(1.5, n), rep(1.5, n) )
  P$sense <- c("=", rep("<=", n), rep(">=", n))
  
  w_eqw <- rep(1/n, n)
  PS <- PolytopeSampler$new( polytope = P )
  PS$setCtrl( n_sim = 10^5,
              algo = "hitnrun",
              interior_point = w_eqw )
  PS_prime <- PS$copy()
  PS_prime$transform()
  PS_prime$sample()  
  
  # trans <- PS_prime$transformation
  # tmat_inv <- trans$tmat_inv
  # z <- trans$z
  # P_prime <- P$transform( tmat_inv = tmat_inv, z = z, idx_eq = 1 )
  # P_prime  
  # P
  # 
  # samples <- PS_prime$samples$S
  # weights <- t(apply( samples, 1, function(y) { tmat_inv %*% (c(0, y) + z) } ) )
  # apply( weights, 1, sum )
  
  PS_new <- PS_prime$copy()
  PS_new$transformBack()
  weights <- PS_new$samples$S  
  
  
  theta_ratio <- apply( weights, 1, function(w) { (w %*% z2) / (w %*% z1) } )
  qvec_ratio <- quantile( theta_ratio, seq(from = 0, to = 1, length.out = nq) )
  colors_ratio <- theta_ratio * NA
  for ( i in length(qvec_ratio):2 ) {
    colors_ratio[ which( theta_ratio <= qvec_ratio[i] ) ] <- color_vec_ratio[i]
  }
  
  
  S <- Simplex$new( d = n )
  weights_simplex <- S$runif( n_sim = n_sim )
  theta_ratio_simplex <- apply( weights_simplex, 1, function(w) { (w %*% z2) / (w %*% z1) } )
  qvec_ratio_simplex <- quantile( theta_ratio_simplex, seq(from = 0, to = 1, length.out = nq) )
  colors_ratio_simplex <- theta_ratio_simplex * NA
  for ( i in length(qvec_ratio_simplex):2 ) {
    colors_ratio_simplex[ which( theta_ratio_simplex <= qvec_ratio_simplex[i] ) ] <- color_vec_ratio[i]
  }
  
  plot3d( x = weights[ ,1], y = weights[ ,2], z = weights[ ,3], col = colors_ratio, size = 2 )
  points3d( x = weights_simplex[ ,1], y = weights_simplex[ ,2], z = weights_simplex[ ,3], col = colors_ratio_simplex, size = 2 )
  

  
  
  
  PS$setCtrl( n_sim = 10^5,
              algo = "volesti",
              interior_point = w_eqw )
  PS_prime <- PS$copy()
  PS_prime$transform()
  # debugonce( PS_prime$sample )
  PS_prime$sample()  
  
  
  
  
  