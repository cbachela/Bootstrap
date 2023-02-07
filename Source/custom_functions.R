  
  
  ############################################################################
  ### CUSTOM FUNCTIONS, BOOTSTRAP
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     30.03.2020
  # First version:    23.03.2020
  # --------------------------------------------------------------------------
  
  # meanBoot
  # permutations
  # lpNorm
  # comb
  # combMpfr
  # combMpfr2
  # gridDraws
  # draws2prob
  # diffentropyDirichlet
  # skewFUN
  # skew.wt
  # M3BFUN
  # M3BBFUN
  # kurtFUN
  # kurt.wt
  # M4BFUN
  # M4BBFUN
  # momFUN
  # betaMoments
  # getMk
  # test.getMk
  # varsi
  # varsiFindQuantile
  # lasserre
  # biasCorrection
  # M2exact
  # M2Analytic
  # M3Exact
  # M4Exact
  # mn2Lambda
  # kernelFunction
  
  
  
  
  
  # --------------------------------------------------------------------------
  biasCorrection <- function( N = NULL, w_bar = NULL, m = 2, unweighted = FALSE )
  {
    # see Rimoldini_2013 Weighted skewness and kurtosis unbiased by sample size
    
    if ( isTRUE(unweighted) ) {
      
      if ( is.null(N) ) {
        stop("'N' needs to be provided.")
      }
      if ( m == 2 ) {
        ans <- N / (N - 1)
      } else if ( m == 3 ) {
        ans <- N^2 / ((N - 1) * (N - 2))
      } else if ( m == 4 ) {
        scl4denom <- (N - 1) * (N - 2) * (N - 3)
        scl4a <- N * (N^2 - 2 * N + 3)
        scl4b <- 3 * N * (2 * N - 3)
        ans <- c( scl4a = scl4a, scl4b = scl4b, scl4denom = scl4denom )
      } else {
        stop( "parameter m has to be 2, 3 or 4.")
      }
      
    } else {
      
      if ( is.null(w_bar) ) {
        if ( is.null(N) ) {
          stop("Either 'N' or 'w_bar' needs to be provided.")
        }
        w_bar <- rep(1/N, N)
      }
      
      V1 <- sum(w_bar)
      V2 <- sum(w_bar^2)
      V3 <- sum(w_bar^3)
      V4 <- sum(w_bar^4)
      
      if ( m == 2 ) {
        ans <- V1^2 / (V1^2 - V2)
      } else if ( m == 3 ) {
        ans <- V1^3 / (V1^3 - 3 * V1 * V2 + 2 * V3)
      } else if ( m == 4 ) {
        scl4denom <- ((V1^2 - V2) * (V1^4 - 6*V1^2*V2 + 8*V1*V3 + 3*V2^2 - 6*V4))
        scl4a <- (V1^2 * (V1^4 - 3*V1^2*V2 + 2*V1*V3 + 3*V2^2 - 3*V4))
        scl4b <- (3*V1^2 * (2*V1^2*V2 - 2*V1*V3 - 3*V2^2 + 3*V4))
        ans <- c( scl4a = scl4a, scl4b = scl4b, scl4denom = scl4denom )
      } else {
        stop( "parameter m has to be 2, 3 or 4.")
      }
    }
    
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  M2Exact <- function( z = NULL, 
                       scl2 = NULL, 
                       m2 = NULL, 
                       w_bar = NULL, 
                       exp2norm2 = NULL,
                       M = NULL )
  {
    stopifnot( is.numeric(exp2norm2) )
    if ( is.null(w_bar) ) {
      N <- length(z)
      w_bar <- rep(1/N, N)
    }
    if ( !is.null(M) ) {
      z <- t(M) %*% z
    }
    if ( is.null(m2) ) {
      z_bar <- sum( w_bar * z )
      m2_ml <- sum( w_bar * (z - z_bar)^2 )
      if ( is.null(scl2) ) {
        scl2 <- biasCorrection( w_bar = w_bar, m = 2 )
      }
      m2 <- m2_ml * scl2
    }
    ans <- m2 * ( exp2norm2 - sum(w_bar^2) )
    return( ans )
  }
  # --------------------------------------------------------------------------
  M2Analytic <- function( z = NULL, 
                          w_bar = NULL, 
                          m2_ml = NULL,
                          method = c("classic", "Bayesian") )
  {
    method <- match.arg(method)
    
    if ( is.null(w_bar) ) {
      N <- length(z)
      w_bar <- rep(1/N, N)
    }
    if ( is.null(m2_ml) ) {
      z_bar <- sum( w_bar * z )
      m2_ml <- sum( w_bar * (z - z_bar)^2 )
    }
    
    if ( method == "classic" ) {
      ans <- m2_ml / n
    } else {
      ans <- m2_ml / (n + 1)
    }
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  M3Exact <- function( z = NULL, 
                       scl3 = NULL, 
                       m3 = NULL, 
                       w_bar = NULL, 
                       exp3norm3 = NULL, # E(||w||_3^3)
                       # exp2norm2 = NULL, # E(||w||_2^2)
                       expw2 = NULL, # E(w^2)
                       M = NULL )
  {
    if ( is.null(w_bar) ) {
      N <- length(z)
      w_bar <- rep(1/N, N)
    }
    if ( !is.null(M) ) {
      z <- t(M) %*% z
    }
    if ( is.null(m3) ) {
      z_bar <- sum( w_bar * z )
      m3_ml <- sum( w_bar * (z - z_bar)^3 )
      if ( is.null(scl3) ) {
        scl3 <- biasCorrection( w_bar = w_bar, m = 3 )
      }
      m3 <- m3_ml * scl3
    }
    # ans <- m3 * ( exp3norm3 - 3 * exp2norm2 * sum(w_bar^2) + 2 * sum(w_bar^3) ) # only holds if w_bar = 1/n, for all i
    ans <- m3 * ( exp3norm3 - 3 * sum(w_bar * expw2) + 2 * sum(w_bar^3) ) # only holds if w_bar = 1/n, for all i
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  M4Exact <- function( z = NULL, 
                       scl4 = NULL, 
                       m4 = NULL, 
                       w_bar = NULL, 
                       exp4norm4 = NULL, # E(||w||_4^4)
                       expw3 = NULL, # E(w^3)
                       expw2 = NULL, # E(w^2)
                       M = NULL )
  {
    
    # WRONG !!!
    
    if ( is.null(w_bar) ) {
      N <- length(z)
      w_bar <- rep(1/N, N)
    }
    if ( !is.null(M) ) {
      z <- t(M) %*% z
    }
    if ( is.null(m4) ) {
      z_bar <- sum( w_bar * z )
      m2_ml <- sum( w_bar * (z - z_bar)^2 )
      m4_ml <- sum( w_bar * (z - z_bar)^4 )
      if ( is.null(scl4) ) {
        scl4 <- biasCorrection( w_bar = w_bar, m = 4 )
      }
      m4 <- (m4_ml * scl4["scl4a"] - m2_ml^2 * scl4["scl4b"]) / scl4["scl4denom"]
    }
    # ans <- m4 * ( exp4norm4 - 3 * sum(w_bar^4) - 4 * exp3norm3 * sum(w_bar^2) + 6 * exp2norm2 * sum(w_bar^3) )
    ans <- m4 * ( exp4norm4 - 3 * sum(w_bar^4) - 4 * sum(w_bar * expw3) + 6 * sum(w_bar^2 * expw2) )
    return( ans )
  }
  
  
  
  
  

  # --------------------------------------------------------------------------
  meanBoot <- function(x, idx)
  {
    if ( NCOL(x) > 1 ) {
      ans <- apply(x[idx, ], 2, mean)
    } else {
      ans <- mean(x[idx])
    }
    return( ans )
  }
  
  # FROM RP PACKAGE:
  # --------------------------------------------------------------------------
  permutations <- function(x) 
  {
    if (length(x) == 1) {
      return(x)
    }
    else {
      res <- matrix(nrow = 0, ncol = length(x))
      for (i in seq_along(x)) {
        res <- rbind(res, cbind(x[i], Recall(x[-i])))
      }
      return(res)
    }
  }
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  lpNorm <- function(x, x_bar = rep(0, length(x)), p = 1, b_abs = TRUE ) 
  { 
    z <- x - x_bar
    if ( isTRUE(b_abs) ) {
      z <- abs(z)
    }
    ans <- sum(z^p)^(1/p) 
    return( ans )
  } 
  
  # --------------------------------------------------------------------------
  # comb <- slolz:::comb
  comb <- function(n, k) { factorial(n) / ( factorial(n-k) * factorial(k)) }
  combMpfr <- function(n, k) { Rmpfr::chooseMpfr(n, k) }
  combMpfr2 <- function(n, k) { Rmpfr::factorialMpfr(n) / ( Rmpfr::factorialMpfr(n-k) * Rmpfr::factorialMpfr(k)) }
  
  
  # --------------------------------------------------------------------------
  gridDraws <- function ( parts_mat, cnames = FALSE, mpfr = FALSE ) 
  {
  
    if ( isTRUE(mpfr) ) {
      FUN <- function ( x ) 
      {
        tmp <- table(x)
        num1 <- Rmpfr::factorialMpfr(sum(x))
        if ( length(x) == sum(x) ) {
          num2 <- num1
        } else {
          num2 <- Rmpfr::factorialMpfr(length(x))
        }
        n_paths_to_level_set <- num1 / prod( Rmpfr::factorialMpfr(x) )
        n_points_in_level_set <- num2 / prod( Rmpfr::factorialMpfr(tmp) )
        ans <- c(n_paths_to_level_set,
                 n_points_in_level_set)
        return( ans )
      }
    } else {
      FUN <- function ( x ) 
      {
        tmp <- table(x)
        num1 <- factorial(sum(x))
        if ( length(x) == sum(x) ) {
          num2 <- num1
        } else {
          num2 <- factorial(length(x))
        }
        n_paths_to_level_set <- num1 / prod(factorial(x))
        n_points_in_level_set <- num2 / prod(factorial(tmp))
        ans <- c(n_paths_to_level_set,
                 n_points_in_level_set)
        return( ans )
      }
    }
    ans <- apply( X = parts_mat, MARGIN = 2, FUN = FUN )
    if ( isTRUE(mpfr) ) {
      ans <- do.call( cbind, ans )
    }
    rownames(ans) <- c("# paths to the same point", 
                       "# points in the same level set")
    if ( isTRUE(cnames) ) {
      colnames(ans) <- apply( X = parts_mat, 
                              MARGIN = 2, 
                              FUN = function(X) { paste(X, collapse = ".") } )
    }
    
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  draws2prob <- function ( draws ) 
  {
    ans <- draws[1, ] / as.numeric(draws[1, ] %*% draws[2, ])
    return( ans )
  }
  
  
  
  # --------------------------------------------------------------------------
  # Differential entropy of Dirichlet symmetric random variable 
  # --------------------------------------------------------------------------
  diffentropyDirichlet <- function( alpha )
  {
    if ( all(alpha == 1) ) {
      H <- log( 1 / gamma( n ) )
    } else {
      n <- length(alpha)
      a0 <- sum(alpha)
      B <- prod( gamma(alpha) ) / gamma( a0 )
      H <- log(B) + (a0 - n) * digamma(a0) - sum( (alpha - 1) * digamma(alpha) )
    }
    return( H )
  }
  
  
  
  # --------------------------------------------------------------------------
  # Skewness 
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  skewFUN <- function(x) { n <- length(x); sum( (x - mean(x))^3 ) / (n * (var(x) * (n-1) / n)^(3/2) ) }
  
  # --------------------------------------------------------------------------
  skew.wt <- function(x, wt) 
  {
    mu_wt <- sum( wt * x )
    sigma_wt <- cov.wt( x = matrix(x, ncol = 1), 
                        wt = wt, 
                        cor = FALSE, 
                        center = TRUE, 
                        method = "ML" )$cov
    ans <- sum( wt * (x - mu_wt)^3 ) / sigma_wt^(3/2)
    return( as.numeric(ans) )
  }
  
  # --------------------------------------------------------------------------
  M3BBFUN <- function( x, M2 = NULL, wghts = NULL, scl = TRUE )
  { 
    # Skewness of exact Bayesian bootstrap distribution of mean
    n <- length(x)
    if ( is.null(wghts) ) {
      mu <- mean(x)
    } else {
      mu <- sum( wghts * x )
    }
    if ( !is.null(M2) ) {
      if ( is.null(wghts) ) {
        ans <- (2 / (n * (n + 1) * (n + 2))) * 
                  sum( (x - mu)^3 )
      } else {
        ans <- (2 / ((n + 1) * (n + 2))) * 
                  sum( wghts * (x - mu)^3 ) 
      }
      if ( isTRUE(scl) ) {
        ans <- ans / M2^(3/2)
      }
    } else {
      if ( is.null(wghts) ) {
        sk <- skewFUN( x = x )
      } else {
        sk <- skew.wt( x = x, wt = wghts )
      }
      if ( isTRUE(scl) ) {
        ans <- sk * (2 * sqrt(n + 1)) / (n + 2)
      } else {
        stop("not implemented.")
      }
    }
    return( as.numeric(ans) )
  }
  
  
  # --------------------------------------------------------------------------
  M3BFUN <- function( x, M2 = NULL, wghts = NULL, scl = TRUE )
  { 
    # Skewness of exact classical bootstrap distribution of mean
    n <- length(x)
    if ( is.null(wghts) ) {
      mu <- mean(x)
    } else {
      mu <- sum( wghts * x )
    }
    if ( !is.null(M2) ) {
      if ( is.null(wghts) ) {
        ans <- (1 / M2^(3/2)) * (1 / n^3) * sum( (x - mu)^3 )
      } else {
        ans <- (1 / M2^(3/2)) * (1 / n^2) * sum( wghts * (x - mu)^3 )
      }
    } else {
      if ( is.null(wghts) ) {
        sk <- skewFUN( x = x )
      } else {
        sk <- skew.wt( x = x, wt = wghts )
      }
      if ( isTRUE(scl) ) {
        ans <- sk * n^(3/2) / n^2
      } else {
        if ( is.null(wghts) ) {
          wghts <- rep(1/n, n)
        }
        ans <- sum( wghts * (x - mu)^3 ) / n^2
      }
    }
    return( as.numeric(ans) )
  }
  
  
  # --------------------------------------------------------------------------
  # Kurtosis
  # --------------------------------------------------------------------------
  kurtFUN <- function(x) { n <- length(x); sum( (x - mean(x))^4 ) / (n * (var(x) * (n-1) / n)^(4/2) ) }
  
  # --------------------------------------------------------------------------
  kurt.wt <- function(x, wt) 
  {
    mu_wt <- sum( wt * x )
    sigma_wt <- cov.wt( x = matrix(x, ncol = 1), 
                        wt = wt, 
                        cor = FALSE, 
                        center = TRUE, 
                        method = "ML" )$cov
    ans <- sum( wt * (x - mu_wt)^4 ) / sigma_wt^(4/2)
    return( as.numeric(ans) )
  }
  
  # --------------------------------------------------------------------------
  M4BBFUN <- function( x, M2 = NULL, wghts = NULL, scl = TRUE ) 
  { 
    # Kurtosis of exact Bayesian bootstrap distribution of mean
    n <- length(x)
    if ( is.null(wghts) ) {
      mu <- mean(x)
    } else {
      mu <- sum( wghts * x )
    }
    if ( !is.null(M2) ) {
      if ( is.null(wghts) ) {
        # ans <- (1 / (n * (n + 1) * (n + 2) * (n + 3))) *
        #           ( 6 * sum( (x - mu)^4 ) + 3 * sum( (x - mu)^2 )^2 )
        m2 <- sum( (x - mu)^2 ) / n
        m4 <- sum( (x - mu)^4 ) / n
      } else {
        # ans <- (1 / ((n + 1) * (n + 2) * (n + 3))) *
        #          ( 6 * sum( wghts * (x - mu)^4 ) + 3 * n * ( sum( wghts * (x - mu)^2 ) )^2 )
        
        m2 <- sum( wghts * (x - mu)^2 )
        m4 <- sum( wghts * (x - mu)^4 )
      }
      ans <- ( 6 * m4 + 3 * n * m2^2 ) / ( (n + 1) * (n + 2) * (n + 3) )
      if ( isTRUE(scl) ) {
        ans <- ans / M2^2
      }
    } else {
      if ( is.null(wghts) ) {
        ks <- momFUN( x = x, k = 4, scl = scl )
      } else {
        ks <- momFUN.wt( x = x, k = 4, wt = wghts, scl = scl )
      }
      if ( isTRUE(scl) ) {
        ans <- (2 * ks + n) * (3 * (n + 1)) / ((n + 2) * (n + 3))
      } else {
        stop("not implemented.")
      }
    }
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  M4BFUN <- function( x, M2 = NULL, wghts = NULL, scl = TRUE ) 
  { 
    # Kurtosis of exact classical bootstrap distribution of mean
    n <- length(x)
    if ( is.null(wghts) ) {
      mu <- mean(x)
    } else {
      mu <- sum( wghts * x )
    }
    if ( !is.null(M2) ) {
      if ( !isTRUE(scl) ) {
        stop("incompatible parameter specifications: scl cannot be FALSE if M2 is provided.")
      }
      if ( is.null(wghts) ) {
        ans <- (1 / M2^2) * (1 / n^4) * sum( (x - mu)^4 ) + 3 - 3 / n
      } else {
        ans <- (1 / M2^2) * (1 / n^3) * sum( wghts * (x - mu)^4 ) + 3 - 3 / n
      }
    } else {
      if ( is.null(wghts) ) {
        ks <- momFUN( x = x, k = 4, scl = scl )
      } else {
        ks <- momFUN.wt( x = x, k = 4, wt = wghts, scl = scl )
      }
      if ( isTRUE(scl) ) {
        ans <- ks / n + 3 - 3 / n
      } else {
        if ( is.null(wghts) ) {
          wghts <- rep(1/n, n)
        }
        ans <- 3 * sum( wghts * (x - mu)^2 )^2 / n^2 + 
                ( sum( wghts * (x - mu)^4 ) - 3 * sum( wghts * (x - mu)^2 )^2 ) / n^3
      }
    }
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  momFUN <- function(x, k, scl = TRUE) 
  { 
    n <- length(x)
    ans <- sum( (x - mean(x))^k ) / n 
    if ( isTRUE(scl) ) {
      ans <- ans / (var(x) * (n-1) / n)^(k/2)
    }
    return( ans )
  }
  # --------------------------------------------------------------------------
  momFUN.wt <- function(x, k, wt, scl = TRUE) 
  {
    mu_wt <- sum( wt * x )
    ans <- sum( wt * (x - mu_wt)^4 )
    if ( isTRUE(scl) ) {
      sigma_wt <- cov.wt( x = matrix(x, ncol = 1), 
                          wt = wt, 
                          cor = FALSE, 
                          center = TRUE, 
                          method = "ML" )$cov
      ans <- ans / sigma_wt^(k/2)
    }
    return( as.numeric(ans) )
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  MkDirichlet <- function( x, M2 = NULL, alpha = NULL ) 
  { 
    # k-th moment of exact Bayesian bootstrap distribution of mean
    
  }
  
  # --------------------------------------------------------------------------
  MkNorm <- function( x, M2 = NULL, x0 = NULL, e2n = NULL ) 
  { 
    # k-th moment of exact bootstrap distribution of mean 
    
    
  }
  
  
  
  
  # --------------------------------------------------------------------------
  betaMoments <- function( a, a0, m )
  {
    mvec <- 0:(m-1)
    # prod( a + mvec ) / prod( a + b + mvec )
    prod( (a + mvec) / (a0 + mvec) )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  getMk <- function( z, k, M2 = NULL, wghts = NULL )
  {
    
    # // Algorithm 3 in CalesChalkisEmiris2021
    
    if ( is.null(wghts) ) {
      wghts <- rep(1/length(z), length(z))
    } 
    z_bar <- sum( wghts * z )
    
    # 1. Compute Av as in Lemma 14: Av_i = R_i - M1
    Av <- z - z_bar
    
    if ( is.null(M2) ) {
      # 2. Compute M_2 by Theorem 15
      M2 <- sum( wghts * Av^2 ) / (n + 1)
    }
    
    # 3. Compute partitions
    parts_mat <- partitions::parts( k )
    
    # 4. Set S = 0
    S <- 0
    
    # 5. For each partition:
    for ( j in 1:ncol(parts_mat) ) {
      
      # (a) decompose the partitions in its d non-zero elements (l_i)_{i=1}^d with multiplicities (p_i)_{i=1)^d
      l_tmp <- sort( parts_mat[parts_mat[ ,j] > 0, j] )
      l <- unique( l_tmp )
      p <- unlist( lapply( l, FUN = function(i) { length(which(l_tmp == i)) } ) )
      
      # (b) a = \prod_{i=1}^d p_i! l_i^{p_i}
      a <- prod( factorial(p) * l^p )
      
      # (c) b = \prod_{i=1}^d ( \sum_{j=1}^n Av_j^{l_i} )^{p_i}
      b <- 1
      for ( i in 1:length(l) ) {
        # b <- b * sum( wghts^2 ) * ( sum( Av^l[i] ) )^p[i]
        # b <- b * sum( wghts^p[i] ) * ( sum( Av^l[i] ) )^p[i]
        # b <- b * sum( wghts^(k/2) ) * ( sum( Av^l[i] ) )^p[i]
        b <- b * ( sum( Av^l[i] ) )^p[i]
      }
      
      # (d) S = S + a / b
      # S <- S + a / b   # error in paper
      S <- S + b / a
      
    }
    
    # 6. Set M_k = S / ( sqrt(M2)^k \choose(n - 1 + k, k) )
    M_k <-  1 / M2^(k/2) * S / choose(n = n - 1 + k, k)
    attr(M_k, "M2") <- M2
    attr(M_k, "S") <- S
    attr(M_k, "B") <- S / choose(n = n - 1 + k, k)
    
    return( M_k )
  }
  
  
  
  # --------------------------------------------------------------------------
  test.getMk <- function()
  {
    
    n <- 10^2
    n_sim <- 10^3
    n_boot <- 50
    set.seed(1111)
    df <- 5
    x <- rt( n = n, df = df )
    z <- x^2
    
    wghts <- rep(1/n, n)
    M2 <- sum( wghts * (z - sum(wghts * z))^2 ) / (n + 1)
    
    lMk <- lapply( 3:5, FUN = function(k) { getMk( z = z, k = k ) } )
    Mk <- unlist( lMk )
    S <- unlist( lapply( lMk, FUN = function(x) { attr(x, "S") } ) )
    B <- unlist( lapply( lMk, FUN = function(x) { attr(x, "B") } ) )
    
    cbind( Mk, 1 / M2^((3:5)/2) * B )  
    
    

    # Run bootstrap
    boot_mat <- matrix( NA, nrow = n_sim, ncol = n_boot )
    bb_mat <- boot_mat
    for ( j in 1:n_boot ) {
      Boot <- boot::boot( data = z,
                          statistic = meanBoot,
                          weights = wghts,
                          R = n_sim )
      boot_mat[ ,j] <- Boot$t
      lP[[j]] <- rdirichlet( n = n_sim, alpha = wghts * n )
      bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
    }
    
    
    # Central moments (Angelova)
    mu <- unlist( lapply( 1:5, FUN = function(i) { sum( (z - mean(z))^i ) / n } ) )
    m1 <- mean(z)
    m2 <- sum( (z - m1)^2 ) / n^2
    m3 <- sum( (z - m1)^3 ) / n^3
    m4 <- 3 * mu[2]^2 / n^2 + (mu[4] - 3 * mu[2]^2) / n^3
    m5 <- 10 * mu[3] * mu[2] / n^3 + (mu[5] - 10 * mu[3] * mu[2]) / n^4
    mom_angelova <- c(m3, m4, m5)
    
    
    
    cbind( mom_angelova, B )
    
    sum( (z - m1)^3 ) / n^3
    sum( (z - m1)^3 ) * 2 / (n * (n + 1) * (n + 2))
    
    mean( apply( boot_mat, 2, function(x) { sum( (x - mean(x))^3 ) / length(x) } ) )
    mean( apply( bb_mat, 2, function(x) { sum( (x - mean(x))^3 ) / length(x) } ) )
    
    
    
    S / B
    choose(n = n - 1 + 3:5, 3:5)
    
    S / mom_angelova
    choose(n = n - 1 + 3:5, 3:5)

    
    
    
    mom_angelova / choose(n = n - 1 + 3:5, 3:5)
    
    
    
  }
  
  
  
  # --------------------------------------------------------------------------
  #' Varsi's algorithm
  #' @export
  # --------------------------------------------------------------------------
  varsi <- function( mu, b )
  {
    z <- mu - b
    idx_neg <- which(z < 0)
    idx_nonneg <- which(z >= 0)
    J <- length(idx_neg)
    K <- length(idx_nonneg)
    x <- 0
    y <- 0
    if ( J > 0 ) {
      x <- z[idx_neg]
    }
    if ( K > 0 ) {
      y <- z[idx_nonneg]
    }
    a <- c(1, rep(0, K))
    if ( J > 0 ) {
      A <- matrix( NA, nrow = J, ncol = K+1 )
      for ( j in 1:J ) {
        for ( k in 1:K ) {
          a[k+1] <- (y[k] * a[k+1] - x[j] * a[k]) / (y[k] - x[j])
        }
        A[j, ] <- a
      }
    } else {
      A <- a
    }
    
    return( A )
  }
  
  # varsi <- function( mu, b )
  # {
  #   z <- mu - b
  #   idx_neg <- which(z < 0)
  #   idx_nonneg <- which(z >= 0)
  #   J <- length(idx_neg)
  #   K <- length(idx_nonneg)
  #   if ( J > 0 ) {
  #     x <- z[idx_neg]
  #   }
  #   if ( K > 0 ) {
  #     y <- z[idx_nonneg]
  #   }
  #   a <- c(1, rep(0, K))
  #   if ( J > 0 ) {
  #     for ( j in 1:J ) {
  #       for ( k in 1:K ) {
  #         a[k+1] <- (y[k] * a[k+1] - x[j] * a[k]) / (y[k] - x[j])
  #       }
  #     }
  #   }
  # 
  #   return( tail(a, 1) )
  # }
  
  
  # --------------------------------------------------------------------------
  varsiFindQuantile <- function( z, th, tol = 1e-12 )
  {
    FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) }
    z0 <- quantile(z, th)
    a <- min(z)
    b <- max(z)
    i <- 1
    err <- 1
    tic <- Sys.time()
    while( err > tol ) {
      p <- FUN( z0 = z0 )
      err <- abs(p - th)
      if ( p > th ) {
        b <- z0
      } else {
        a <- z0
      }
      z0 <- (a + b) / 2
      i <- i + 1
    }
    toc <- Sys.time() - tic
    attr(z0, "i") <- i
    attr(z0, "tictoc") <- toc
    return( z0 )   
  }
  
  # --------------------------------------------------------------------------
  #' Lasserre's analytic solution
  #' @description Implements Theorem 2.2 in Lasserre (2015)
  #' @references J. B. Lasserre (2015). Volume of slices and section of the simplex in closed form
  #' @export
  # --------------------------------------------------------------------------
  lasserre <- function( a, th )
  {
    
    N <- length(a)
    a <- c(0, a)
    ans <- 0
    for ( i in 1:length(a) ) {
      num <- max(0, th - a[i])^N
      denom <- prod( a[-i] - a[i] )
      ans <- ans + num / denom  
    }
    # ans <- ans  / factorial(N)
    return( ans )    
  }
  
  
  test.lasserre <- function()
  {
    require(RP)
    set.seed(1234)
    n <- 100
    x <- rnorm(n)
    b <- quantile(x, 0.6)
    
    tail( as.numeric( varsi( mu = x, b = b ) ), 1 )
    lasserre( a = x, th = b ) 
    lasserre( a = x, th = b ) * factorial(n)
    
    
    
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Function finds lambda such that the variance of the Bayesian bootstrap
  # based on Dir(\alpha \lambda) corresponds to the variance of the m out of n bootstrap.
  # --------------------------------------------------------------------------
  mn2Lambda <- function(n, m) { (1-m)/(1-n) + (m-1)/(n*(1-n)) }
  
  
  
  
  # --------------------------------------------------------------------------
  kernel.gaussian <- function(x) { 1 / sqrt(2 * pi) * exp( -x^2 / 2 ) }
  # --------------------------------------------------------------------------
  kernel.epanechnikov <- function(x) { 3 / 4 * (1 - x^2) * (abs(x) < 1) }
  # --------------------------------------------------------------------------
  kernel.boxcar <- function(x) { 1 / 2 * (abs(x) < 1) }
  # --------------------------------------------------------------------------
  kernel.tricube <- function(x) { 70 / 81 * (1 - abs(x)^3) * (abs(x) < 1) }
  # --------------------------------------------------------------------------
  kernelFunction <- function(x, method = c("gaussian", "epanechnikov", "boxcar", "tricube") )
  {
    method <- match.arg(method)
    FUN <- match.fun( paste0("kernel.", method) )
    y <- FUN(x = x )
    
    return( y )
  }
  
  
  
  
  
  
  
  
  
  