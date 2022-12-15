  
  
  ############################################################################
  ### REGRESSION BOOTSTRAP - RESIDUALS RESAMPLING
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     23.06.2022
  # First version:    23.06.2022
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
  
  # TD <- trainingData( Y_train = cbind(X_bm, X),
  #                     X_train = cbind(FF3, FF5) )
  
  ###
  TD <- trainingData( Y_train = head( X_bm, 10 ),
                      X_train = X[ ,1:3] )
  TD
  ###
  
  
  
  
  
  # Synthetic data 
  n <- 10^2
  set.seed(1111)
  a <- 0.01
  b <- 0.2
  # Normal distribution
  y <- scale( rnorm( n, 0, 1 ), TRUE, TRUE )
  r <- scale( rnorm( n, 0, 1 ), TRUE, TRUE )
  # Student t
  y <- scale( rt( n, df = 4 ), TRUE, TRUE )
  r <- scale( rt( n, df = 4 ), TRUE, TRUE )
 
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
  b_hat <- coeff[2:nrow(coeff), 1]
  resid <- residuals(reg$reg)
  summary(reg$reg)  
  

  # debugonce( betaHat )
  betaHat( X = X_train, y = Y_train, wghts = NULL, intercept = TRUE, type = "ols" )
  coeff
  
  
  
  # --------------------------------------------------------------------------
  # Bootstrap (classical) of residuals, X fixed
  #
  # (assumes that the functional form of the regression model fit to the data is correct 
  # and that the errors are identically distributed
  # --------------------------------------------------------------------------
 
  
  n_sim <- 10^3
  n <- nrow(TD$Y_train)
  gridnodes <- matrix( 0, nrow = n_sim, ncol = n )
  ab_boot <- matrix( 0, nrow = n_sim, ncol = ncol(X_train)+1, 
                     dimnames = list(NULL, c("a", paste0("b", 1:ncol(X_train)))) )
  ab_boot2 <- ab_boot
  ab_analyt <- ab_boot
  
  for ( i in 1:n_sim ) {
    
    idx <- sample( x = 1:n, size = n, replace = TRUE )
    resid_boot <- resid[idx]
    y_boot <- timeSeries( a_hat + X_train %*% b_hat + resid_boot, rownames(X_train) )
    reg_tmp <- regression( Y_train = y_boot,
                           X_train = X_train,
                           type = "ols" )
    ab_boot[i, ] <- reg_tmp$coeffmat[ ,1]
    
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    wghts <- gridnodes[i, ]
    y_boot2 <- timeSeries( a_hat + X_train %*% b_hat + wghts * n * resid, rownames(X_train) )
    reg_tmp2 <- regression( Y_train = y_boot2,
                            X_train = X_train,
                            type = "ols" )
    ab_boot2[i, ] <- reg_tmp2$coeffmat[ ,1]
    
    # debugonce( alphaBetaHatBoot )
    ab_analyt[i, ] <- alphaBetaHatBoot( x = X_train,
                                        y = Y_train,
                                        w = wghts,
                                        a_hat = a_hat,
                                        b_hat = b_hat,
                                        resid = resid )
    
  }
  
 
  a_mat <- cbind( ab_boot[ ,"a"], ab_boot2[ ,"a"], ab_analyt[ ,"a"] )
  b1_mat <- cbind( ab_boot[ ,"b1"], ab_boot2[ ,"b1"], ab_analyt[ ,"b1"] )
  # b2_mat <- cbind( ab_boot[ ,"b2"], ab_boot2[ ,"b2"], ab_analyt[ ,"b2"] )
  
  
  coeff[1, ]
  apply( a_mat, 2, mean ); apply( a_mat, 2, sd )
  barplot( abs( coeff[1, 1] - apply( a_mat, 2, mean ) ) )
  barplot( abs( coeff[1, 2] - apply( a_mat, 2, sd ) ) )
  
  coeff[2, ]
  apply( b1_mat, 2, mean ); apply( b1_mat, 2, sd )
  barplot( abs( coeff[2, 1] - apply( b1_mat, 2, mean ) ) )
  barplot( abs( coeff[2, 2] - apply( b1_mat, 2, sd ) ) )
  
  
  
  ldens_a <- apply( a_mat, 2, density )
  ldens_b <- apply( b1_mat, 2, density )
  
  slolz:::plot.ldensity( ldens_a, fillin = FALSE, colors = 1:ncol(a_mat) )
  abline( v = a_hat )
  slolz:::plot.ldensity( ldens_b, fillin = FALSE, colors = 1:ncol(b1_mat) )
  abline( v = b_hat[1] )  
  
  
  
  X <- X_train
  A <- solve(t(X) %*% X) %*% t(X)
  a <- A[1, ]
  z <- a * resid * n
  sqrt( var(z) * ( 2/n - 1/n^2 - 1/n ) )
  apply( b1_mat, 2, sd )
  coeff
  
  apply( b1_mat, 2, mean )
  b_hat + mean(z)
  
  
  
  
  
  
  
  
  
  
  ######################
  
  
  
  n_sim <- 10^3
  n <- nrow(TD$Y_train)
  X <- TD$X_train
  y <- TD$Y_train
  gridnodes <- matrix( 0, nrow = n_sim, ncol = n )
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  b_mat <- matrix( NA, nrow = n_sim, ncol = 4 )
  
  for ( i in 1:nrow(samples) ) {
    
    idx <- sample( x = 1:n, size = n, replace = TRUE )
    resid_boot <- resid[idx]
    y_star <- coeff[1, 1] + coeff[2, 1] * X + resid_boot
    
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl / n
    wghts <- gridnodes[i, ]
    
    # debugonce( betaHat )
    b_mat[i, 1] <- betaHat( X = X, y = y_star, type = "ols" )[1]
    b_mat[i, 2] <- betaHat( X = X, y = y, type = "wls", wghts = wghts )[1]
    b_mat[i, -c(1:2)] <- betaHat( X = X, y = y, type = "wrr", wghts = wghts )
    
  }  
  
  ldens <- apply( b_mat, 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE, col = 1:ncol(b_mat) ) 

  
  boxplot( b_mat )
  head(b_mat)
  
  
  
  
  
  
  # Bayesian bootstrap
  
  n_sim <- 10^3
  n <- nrow(TD$Y_train)
  X <- TD$X_train
  y <- TD$Y_train
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  # ...
  