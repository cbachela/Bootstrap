  
  
  ############################################################################
  ### REGRESSION BOOTSTRAP - HELPER FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.07.2022
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
  
  
  
  
  # --------------------------------------------------------------------------
  # Functions
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  alphaHat <- function( y_bar, b_hat, x_bar ) { y_bar - b_hat %*% x_bar }
  
  
  # --------------------------------------------------------------------------
  betaHat <- function( X, y, 
                       wghts = NULL, 
                       intercept = FALSE,
                       type = c("ols", "wls", "wrr") ) 
  { 
    type <- match.arg(type)
    
    if ( type == "ols" ) { # ordinary least squares
     
      X <- scale( X, TRUE, FALSE )
      y_bar <- mean(y)
      y <- y - y_bar
      # ans <- solve(t(X) %*% X) %*% t(X) %*% y
      b_denom <- t(X) %*% X
      b_num <- t(X) %*% y
      ans <- solve(b_denom) %*% b_num
      if ( isTRUE(intercept) ) {
        a <- alphaHat( y_bar = y_bar,
                       x_bar = attr(X, "scaled:center"), 
                       b_hat = ans )
        ans <- c(a, ans)
      }
      attr(ans, "num") <- b_num
      attr(ans, "denom") <- b_denom

   } else if ( type == "wls" ) { # weighted least squares
      
      # y_bar <- sum( wghts * y )
      y_bar <- mean(y)
      y <- y - y_bar
      # x_bar <- apply( X, 2, function(x) { sum( wghts * x ) } )
      x_bar <- apply( X, 2, mean )
      X <- sweep( x = X, MARGIN = 2, STATS = x_bar, FUN = "-" )
      W <- diag(wghts)
      b_denom <- t(X) %*% W %*% X
      b_num <- t(X) %*% W %*% y
      ans <- solve(b_denom) %*% b_num
      attr(ans, "num") <- b_num
      attr(ans, "denom") <- b_denom
    
    } else if ( type == "wrr" ) { # weighted residuals resampling
      
      X_raw <- X
      x_bar <- apply( X_raw, 2, mean )
      X <- sweep( x = X, MARGIN = 2, STATS = x_bar, FUN = "-" )
      A <- solve( t(X) %*% X ) %*% t(X)
      n <- nrow(X)
      
      beta_ols <- betaHat( X = X_raw, y = y, type = "ols" )[ ,1]
      alpha_ols <- alphaHat( y_bar = mean(y), x_bar = x_bar, b_hat = beta_ols )
      y_ols <- as.numeric( alpha_ols ) + as.numeric( X_raw %*% beta_ols )
      resid <- y - y_ols
      A_tmp <- A * matrix(resid, nrow = nrow(A), ncol = ncol(A), byrow = TRUE) * n
      beta_hat <- A %*% y_ols + A_tmp %*% wghts
      
      beta_wls <- betaHat( X = X_raw, y = y, type = "wls", wghts = wghts )[ ,1]
      alpha_wls <- alphaHat( y_bar = mean(y), x_bar = x_bar, b_hat = beta_wls )
      y_wls <- as.numeric( alpha_wls ) + as.numeric( X_raw %*% beta_wls )
      resid <- y - y_wls
      A_tmp <- A * matrix(resid, nrow = nrow(A), ncol = ncol(A), byrow = TRUE) * n
      beta_hat_2 <- A %*% y_wls + A_tmp %*% wghts
  
      ans <- c( beta_hat, beta_hat_2 )
      
    }
    return( ans )
  }

  
  
  # --------------------------------------------------------------------------
  betaHatUnivariate <- function( x, y, wghts = NULL, type = c("ols", "wls") )
  {
    type <- match.arg(type)
    if ( type == "ols" ) {
      x_bar <- mean(x)
      y_bar <- mean(y)
      b_num <- sum( (y - y_bar) * (x - x_bar) )
      b_denom <- sum( (x - x_bar)^2 )
    } else {
      x_bar <- sum( wghts * x )
      y_bar <- sum( wghts * y )
      b_num <- sum( wghts * (y - y_bar) * (x - x_bar) )
      b_denom <- sum( wghts * (x - x_bar)^2 )
    }
    ans <- b_num / b_denom
    attr(ans, "num") <- b_num
    attr(ans, "denom") <- b_denom
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  alphaBetaHatBoot <- function( x, y, w, a_hat, b_hat, resid )
  {
    # Simple linear regression, i.e., assumes x to be univariate
    n <- NROW(y)
    C <- a_hat + b_hat * x
    C_bar <- mean(C)
    y_boot <- a_hat + b_hat * x + w * n * resid
    y_boot_bar <- C_bar + sum( w * resid )
    x_bar <- mean(x)
    beta_hat_boot <- betaHatUnivariate( x = x, y = y_boot, "ols" )
    # beta_hat_boot <- betaHat( x = x, y = y_boot, x_bar = x_bar, y_bar = y_boot_bar ) 
    alpha_hat_boot <- alphaHat( y_bar = y_boot_bar, 
                                x_bar = x_bar,
                                b_hat = beta_hat_boot )
    
    ans <- c( alpha = alpha_hat_boot,
              beta = beta_hat_boot )
    return( ans )
  }
  
  
  






























