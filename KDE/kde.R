  
  
  ############################################################################
  ### WEIGHTED KERNEL DENSITY ESTIMATION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     27.05.2021
  # First version:    27.05.2020
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(rgl)
  require(boot)
  require(volesti)
  require(slolz)
  require(RP)
  wd <- "H:/R/notingit/Bootstrap/"
  source( paste0(wd, "Source/custom_functions.R") )
     
  
  
  
  
  # --------------------------------------------------------------------------
  # Specifications
  # --------------------------------------------------------------------------
  
  
  n <- 3
  n_sim <- 10^2
  set.seed(1111)
  x <- rnorm(n)
  z <- (x - mean(x))^2
  # names(z) <- aviationAlphabet()[1:n]
  names(z) <- paste0("X", 1:n)
  
  
  
  alpha <- rep(1, n)
  samples <- rdirichlet(n = n_sim, alpha = alpha)
  
  
  require(DAARC)
  
  Obj <- loadSensor( sensor_name = "movavg" )
  X_train <- Obj$getSignal()
  Y_train <- tail( Obj$data$X_bm, 100 )
  n_lag <- 0
  data <- trainingData(Y_train = Y_train,
                       X_train = X_train,
                       n_lag = n_lag)
  
  n <- nrow(data$Y_train)
  # alpha <- rep(1, n)
  # alpha <- rep(1/n, n)
  alpha <- 1:n / sum(1:n) * n
  samples <- rdirichlet(n = n_sim, alpha = alpha)
  ldens <- apply( samples, 1, function(w) { densFUN(data$Y_train, weights = w) } )
  
  slolz:::plot.ldensity( ldens, fillin = FALSE )
  lines(densFUN(data$Y_train), col = 4, lwd = 5 )
  
  
  DS <- dirichletSampling( Y_train = Y_train,
                           X_train = X_train,
                           sclfct = 1,
                           weights_mat = timeSeries(t(samples), time(Y_train)) )  
  ldens_ds <- lapply( DS, FUN = densFUN )
  
  slolz:::plot.ldensity( ldens_ds, fillin = FALSE )
  lines(densFUN(data$Y_train), col = 4, lwd = 5 )
  


  
  stats_fun <- function( p, y )
  { 
    ans <- t(p) %*% y
    return( ans )
  }
  
  
  X_eval <- seq(from = min(data$X_train), 
                  to = max(data$X_train), 
                  length.out = min(10, length(unique(data$X_train))))
  
  # Call function to obtain weighting based on distance to X_eval
  weights_mat <- weightsFun( data, 
                             x_eval = X_eval,
                             method = "kernel" )
  
  
  X_eval
  headleft(weights_mat)
  plot(as.timeSeries(weights_mat))
  
  
  # debugonce( kernel.gaussian )
  n <- nrow(data$Y_train)
  y <- scale(data$Y_train, FALSE, FALSE)
  h <- sd(y)
  kern_y <- kernelFunction( x = y, method = "gaussian", h = h )
  # kern_y <- kernel.gaussian( x = data$Y_train, h = 0.1 )
  plot(kern_y)
  plot( x = y, y = kern_y )
  
  
  dens1 <- density(data$Y_train)
  # dens1 <- densFUN(data$Y_train)
  x <- dens1$x
  bw <- dens1$bw
  
  sum(kern_y) / (h * n)
  
  x0 <- mean(y)
  kernelFunction( x = x0, method = "gaussian", h = h )
  sum( kernelFunction( x = y - x0, method = "gaussian", h = h ) )/ (n*h)
  
  
  tmp <- lapply( dens1$x, function(x0) { sum( kernelFunction( x = y - x0, method = "gaussian", h = bw ) ) / (n*bw) } )
  plot( x = dens1$x, y = unlist(tmp) )
  lines(dens1, col = 2)
  
  
  
  
  
  
  
  
  
  
  if ( is.null(sclfct) ) {
    sclfct <- nrow(data$Y_train)
  }
  if ( isTRUE(scl_by_entropy) ) {
    sclfct <- apply( weights_mat, 2, entropy, exponential = TRUE )
  }
  if ( length(sclfct) == 1 ) {
    sclfct <- rep(sclfct, ncol(weights_mat))
  }
  
  y <- data$Y_train
  if ( ncol(y) == 1 ) {
    y <- as.numeric(y) # for efficiency
  }
  lP <- list()
  lStats <- list()
  for ( j in 1:ncol(weights_mat) ) {
    lP[[j]] <-  rdirichlet(n = n_sim, 
                           alpha = weights_mat[ ,j] * sclfct[j] )
    ###
    if ( sclfct_lp > 1 ) {
      lP[[j]] <- t( apply(lP[[j]], 1, function(x) { x^sclfct_lp / sum(x^sclfct_lp ) }) )
    }
    ###
    lStats[[j]] <- apply( lP[[j]], 
                          1, 
                          stats_fun, 
                          y = y )
  }
  if ( isTRUE(correct_bias) ) {
    mu <- apply( weights_mat, 2, stats_fun, y = y )
    lStats <- lapply( 1:length(mu),
                      FUN = function(i) 
                      { 
                        lStats[[i]] - mean(lStats[[i]]) + mu[i]
                      } )
  }
  
  if ( !is.null(X_eval) ) {
    if ( is.null(dim(X_eval)) ) {
      names(lStats) <- paste0( "X_eval=", round(X_eval, 4) )
    } else {
      names(lStats) <- paste0( "X_eval_", 1:nrow(X_eval) ) 
    }
  } else {
    names(lStats) <- paste0("eval", 1:length(lStats))
  }
  
  attr(lStats, "weights") <- weights_mat
  attr(lStats, "sclfct") <- sclfct
  
  
  
  
  
  

  
  
  
  
    
  
  
  
  
  
  
  
