  
  
  ############################################################################
  ### EXACT STANDARD ERROR
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.05.2021
  # First version:    22.05.2021
  # --------------------------------------------------------------------------
  
  
  # Show that the standard error of the ordinary bootstrap can be computed 
  # exactly.
  
  
  
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
  # Data
  # --------------------------------------------------------------------------
  
  df <- 4
  set.seed(1111)
  x <- rt( n = n, df = df )
  z <- x^2
  
  
  # --------------------------------------------------------------------------
  # Partitions
  # --------------------------------------------------------------------------
  
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  
  # Reorder results by decreasing total probability
  ordering <- rev( order( prob_tot) ) 
  parts_mat <- parts_mat[ ,ordering]
  draws <- draws[ ,ordering]
  prob <- prob[ordering]
  prob_tot <- prob_tot[ordering]
  
  lp_parts <- apply( parts_mat / n, 2, lpNorm, p = 2 )
  
  comb( 2*n-1, n )
  sum( draws[2, ] )
  
  # Case m < n:
  
  parts_mat2 <- partitions::parts( m )
  parts_mat2 <- rbind( parts_mat2, 
                       matrix(0, nrow = n-m, ncol = ncol(parts_mat2)) )
  draws2 <- gridDraws( parts_mat = parts_mat2, cnames = FALSE, mpfr = FALSE )
  prob2 <- draws2prob( draws = draws2 )
  prob_tot2 <- prob2 * draws2[2, ]
  
  # Reorder results by decreasing total probability
  ordering2 <- rev( order( prob_tot2) ) 
  parts_mat2 <- parts_mat2[ ,ordering2]
  draws2 <- draws2[ ,ordering2]
  prob2 <- prob2[ordering2]
  prob_tot2 <- prob_tot2[ordering2]
  
  lp_parts2 <- apply( parts_mat2 / m, 2, lpNorm, p = 2 )
  
  
  
  
  # --------------------------------------------------------------------------
  # Ordinary bootstrap
  # --------------------------------------------------------------------------
  
  boot_mat <- matrix( NA, nrow = B, ncol = n_boot )
  
  tic <- Sys.time()
  for ( j in 1:n_boot ) {
    Boot <- boot::boot( data = z,
                        statistic = meanBoot,
                        R = B )
    boot_mat[ ,j] <- Boot$t
  }
  (toc_boot <- Sys.time() - tic)
  
  
  # Alternatively:
  boot_mat2 <- matrix( NA, nrow = B, ncol = n_boot )
  for ( j in 1:n_boot ) {
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample(1:n, size = m, replace = TRUE)
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl / m
    }
    boot_mat2[ ,j] <- apply( gridnodes, 1, function(p) { sum(p * z) } )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Bayesian bootstrap
  # --------------------------------------------------------------------------
  
  tic <- Sys.time()
  bb_mat <- matrix( NA, nrow = B, ncol = n_boot )
  lP <- list()
  for ( j in 1:n_boot ) {
    lP[[j]] <- rdirichlet( n = B, alpha = rep(1, n) )
    bb_mat[ ,j] <- apply( lP[[j]], 1, function(p) { sum(p * z) } )
  }
  (toc_bb <- Sys.time() - tic)
  
  
  
  
  # --------------------------------------------------------------------------
  # Varsi, i.e., exact BB
  # --------------------------------------------------------------------------
  
  se <- mean( apply(boot_mat, 2, sd) )
  
 
  z_vec <- seq( from = 0, to = mean(z) + se * 10, length.out = 10^5 )
  # FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0) }
  FUN <- function(z0) { tail( as.numeric( varsi( mu = z, b = z0 ) ), 1) }
  tic <- Sys.time()
  p_vec <- try( unlist( lapply( z_vec, FUN ) ) ) 
  if ( inherits(p_vec, "try-error") ) {
    q_vec <- seq( from = 0.1, to = 0.999, length.out = 10^5 )
    z_vec <- quantile(z, q_vec)
    p_vec <- try( unlist( lapply( z_vec, FUN ) ) ) 
  }
  (toc_varsi <- Sys.time() - tic)
  
  # Derivative
  idx <- 3:length(p_vec)
  # d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(q_vec[idx] - q_vec[idx-2])
  d_vec <- (p_vec[idx] - p_vec[idx-2]) / abs(z_vec[idx] - z_vec[idx-2])
  d_vec_std <- d_vec / sum(d_vec)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Standard error
  # --------------------------------------------------------------------------
  
  Vz_unb <- var(z)
  a <- sum( (z - mean(z))^2 )
  Vz_ml <- a / n
  Vz_unb; a / (n - 1)
  Vz_ml
  
  # Ordinary bootstrap
  stderr_boot <- apply( boot_mat, 2, sd )
  stderr_boot2 <- apply( boot_mat2, 2, sd )
  
  # Exact ordinary bootstrap
  stderr_parts <- sqrt( var(z) * ( sum(prob_tot * lp_parts^2) - 1/n ) )
  # stderr_analytic <- sqrt( var(z) * (n-1) / n / n )
  # stderr_analytic <- sqrt( sum( (z - mean(z))^2 ) / n^2 )
  stderr_analytic <- sqrt( Vz_ml / n )
  stderr_parts2 <- sqrt( var(z) * ( sum(prob_tot2 * lp_parts2^2) - 1/n ) )
  
  1 + (sum(prob_tot * lp_parts^2) - 1) / (1 - 1/n)
  1/n # same same !
  
  
 
  
  # Bayesian bootstrap (BB)
  stderr_bb <- apply( bb_mat, 2, sd )
 
  # Exact BB
  stderr_exact_bb <- sqrt( var(z) * ( (2 / (n + 1)) - 1 / n ) )
  stderr_exact_cales <- sqrt( sum((z - mean(z))^2) / (n * (n + 1)) )  # same same
  stderr_exact_cales <- sqrt( var(z) * ((n - 1) / n) / (n + 1) )  # same same
  sqrt( (var(z) * (n - 1)) / (n * (n + 1)) ) # same same
  stderr_exact_bb
  stderr_exact_cales
  sqrt( (sum((z - mean(z))^2) / n) / (n + 1) )  # same same
  Vz_ml * ( 1 + (2/(n+1) - 1) / (1 - 1/n) )
  
  a / (n * (n + 1))
  a * ( 2*n / (n^3 - n) - n / (n^3 - n^2) )
  
  
  1 + (2 / (n + 1) - 1) / (1 - 1/n)
  1 / (n + 1) # same same
  
  
  
   
  # Exact BB by Varsi
  z_eval <- z_vec[-c(1, length(z_vec))]
  sum( d_vec_std * z_eval ); mean(z)
  stderr_exact_bb_varsi <- sqrt( sum( d_vec_std * (z_eval)^2 ) - sum( d_vec_std * z_eval )^2 )
  
  
  
  # Compare
  SE <- setNames( c(mean(stderr_boot), 
                    mean(stderr_boot2),
                    stderr_parts, 
                    stderr_analytic,
                    stderr_parts2,
                    mean(stderr_bb), 
                    stderr_exact_bb,
                    stderr_exact_bb_varsi,
                    sd(z) / sqrt(n)),
                  c("avg_boot",
                    "avg_boot2",
                    "exact_boot", 
                    "exact_boot_analytic",
                    "exact_boot2",
                    "avg_bb", 
                    "exact_bb",
                    "exact_bb_varsi",
                    "normal") )
  SE
  
  
  
  
  # --------------------------------------------------------------------------
  # Skewness
  # --------------------------------------------------------------------------
  
  # Ordinary bootstrap
  M3_boot <- apply( boot_mat, 2, skewFUN )
  M3_boot2 <- apply( boot_mat2, 2, skewFUN )
  
  # Exact ordinary bootstrap
  M3_parts <- M3FUN( x = z, M2 = stderr_parts^2 )  # Do not use, only applies to BB
  M3_parts2 <- M3FUN( x = z, M2 = stderr_parts2^2 ) # Do not use, only applies to BB
  M3_boot_analytic <- sum( (z - mean(z))^3 ) / sum( (z - mean(z))^2 )^(3/2)
  

  
  # Bayesian bootstrap (BB)
  M3_bb <- apply( bb_mat, 2, FUN = skewFUN )
  
  # Exact BB
  M3_exact_bb <- M3FUN( x = z, M2 = stderr_exact_bb^2 )
  M3_exact_bb2 <- M3FUN( x = z, M2 = NULL )
  M3_exact_bb; M3_exact_bb2 # same same
  
  # Compare
  M3 <- setNames( c(mean(M3_boot), 
                    mean(M3_boot2),
                    M3_parts, 
                    M3_parts2,
                    M3_boot_analytic,
                    mean(M3_bb), 
                    M3_exact_bb,
                    0),
                  c("avg_boot",
                    "avg_boot2",
                    "exact_boot", 
                    "exact_boot2",
                    "boot_analytic",
                    "avg_bb", 
                    "exact_bb",
                    "normal") )
  M3
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Kurtosis
  # --------------------------------------------------------------------------
  
  # Ordinary bootstrap
  M4_boot <- apply( boot_mat, 2, kurtFUN )
  M4_boot2 <- apply( boot_mat2, 2, kurtFUN )
  
  # Exact ordinary bootstrap
  M4_parts <- M4FUN( x = z, M2 = stderr_parts^2 )
  M4_parts2 <- M4FUN( x = z, M2 = stderr_parts2^2 )
  mu2 <- sum( (z - mean(z))^2 ) / n
  mu4 <- sum( (z - mean(z))^4 ) / n
  M4_boot_analytic <- (3 * mu2^2 / n^2 + (mu4 - 3 * mu2^2) / n^3) / ( sum( (z - mean(z))^2 / n / n ) )^2
 
  
  # sum( (z - mean(z))^4 ) / (n * (1/n *sum( (z - mean(z))^2 ))^2 )
  # sum( (z - mean(z))^4 ) / sum( (z - mean(z))^2 )^2 * n^2 / n
  # 
  # kurtFUN(z) * n / n^2
  # sum( (z - mean(z))^4 ) / sum( (z - mean(z))^2 )^2
  # mean(M4_boot)
  # 
  # A <- mean(M4_boot) / kurtFUN(z)
  # A
  # n / n^2
  
  
  # Bayesian bootstrap (BB)
  M4_bb <- apply( bb_mat, 2, FUN = kurtFUN )
  
  # Exact BB
  M4_exact_bb <- M4FUN( x = z, M2 = stderr_exact_bb^2 )
  
  
  # Compare
  M4 <- setNames( c(mean(M4_boot), 
                    mean(M4_boot2),
                    M4_parts, 
                    M4_parts2,
                    M4_boot_analytic,
                    mean(M4_bb), 
                    M4_exact_bb,
                    # kurtFUN( rnorm(10^7, mean(z), sd(z) / sqrt(n)))),
                    3),
                  c("avg_boot",
                    "avg_boot2",
                    "exact_boot", 
                    "exact_boot2",
                    "boot_analytic",
                    "avg_bb", 
                    "exact_bb",
                    "normal") )
  M4
  
  
  SE
  M3
  M4
 
  (3 * df - 6) / (df - 4) # Analytic kurtosis of t-dist
  
  
  
  # --------------------------------------------------------------------------
  # Save
  # --------------------------------------------------------------------------
  
  env <- new.env()
  
  env$n <- n
  env$m <- m
  env$B <- B
  env$n_boot <- n_boot
  env$lp_parts <- lp_parts
  env$prob_tot <- prob_tot
  env$SE <- SE
  env$M3 <- M3
  env$M4 <- M4
  
  
  saveRDS( object = env, file = paste0(wd, "waRehouse/exact_std_err_n=", n) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Analyze
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    n_vec <- seq(from = 2, to = 83, by = 1)
    lSE <- list()
    lM3 <- list()
    lM4 <- list()
    e2n <- numeric(length(n_vec))  # expected 2-norm of partitions
    for ( i in seq(along = n_vec) ) {
      env <- readRDS( file = paste0(wd, "waRehouse/exact_std_err_n=", n_vec[i]) )
      lSE[[i]] <- env$SE
      e2n[i] <- sum( env$prob_tot * env$lp_parts^2 )
      lM3[[i]] <- env$M3
      lM4[[i]] <- env$M4
    }
    E2N <- cbind( e2n, 2 / (n_vec + 1), 2 / n_vec )
    E2N <- cbind( E2N, (E2N[ ,2] + E2N[ ,3]) / 2 )
    
    SE <- t( do.call( cbind, lSE ) )
    rownames(SE) <- paste0("n=", n_vec[1:nrow(SE)])
    SE <- cbind( SE, 
                 analytic = unlist( lapply( n_vec, FUN = function(n) { sqrt( var(z[1:n]) * (n-1) / n / n ) } ) ) )
    SE[ ,c("exact_boot", "analytic")]
    E2N
    
    
    plot( e2n )
    lines( 2 / n_vec - 1 / n_vec^2, col = 2 )
    
    
    
    delta <- abs( cbind( E2N[ ,3] - E2N[ ,1], E2N[ ,2] - E2N[ ,1], E2N[ ,4] - E2N[ ,1] ) )
    tail(delta) * 10^3
    plot( as.timeSeries(delta), plot.type = "single" )
    
    
    colors <- c("green", "darkgreen", "red", "darkred", "magenta", "black")
    plot( x = n_vec[1:nrow(SE)], y = SE[ ,1], ylim = range(SE), type = "o" )
    for ( i in 1:ncol(SE) ) {
      lines( x = n_vec[1:nrow(SE)], y = SE[ ,i], col = colors[i] ) 
    }
    
    plot( x = n_vec[1:nrow(SE)], y = (SE[ ,2] - SE[ ,1]), type = "o" )
    lines( x = n_vec[1:nrow(SE)], y = (SE[ ,4] - SE[ ,3]), type = "o", col = 2)
    abline( h = 0, col = "grey" )
    
    
    plot( x = n_vec[1:nrow(SE)], y = (SE[ ,"exact_boot"] - SE[ ,"exact_bb"]), type = "o" )
    
    yrng <- c(0, max( cbind( SE[ ,"normal"] - SE[ ,"exact_boot"], SE[ ,"normal"] - SE[ ,"exact_bb"] ) ) )
    plot( x = n_vec[1:nrow(SE)], y = (SE[ ,"normal"] - SE[ ,"exact_boot"]), type = "o", ylim = yrng )
    lines( x = n_vec[1:nrow(SE)], y = (SE[ ,"normal"] - SE[ ,"exact_bb"]), type = "o", col = 2 )
    abline( h = 0, col = "grey" )  
    
    
    
    
    M3 <-  t( do.call( cbind, lM3 ) )
    M3
    # M3 <- cbind( M3, 
    #              analytic = unlist( lapply( n_vec, FUN = function(n) { skewFUN( z[1:n] ) / n_vec[n]^2 } ) ) )
    # M3[ ,c("avg_boot", "analytic")]
    
    head(M3)
    plot( as.timeSeries(M3[ ,c("avg_boot", "boot_analytic")]), plot.type = "single" )
    plot( as.timeSeries(M3[ ,c("avg_bb", "exact_bb")]), plot.type = "single" )
    
    
    
    M4 <-  t( do.call( cbind, lM4 ) )
    head(M4)
    tail(M4)
    
    plot( as.timeSeries(M4[ ,c("avg_boot", "boot_analytic")]), plot.type = "single" )
    plot( as.timeSeries(M4[ ,c("avg_bb", "exact_bb")]), plot.type = "single" )
    
    
    
    
    
    
    df <- 10
    set.seed(1111)
    x <- rt( n = max(n_vec), df = df )
    z <- x^2
    
    
    tmp <- lapply( 1:nrow(M3), FUN = function(i) { n <- i + 1; sum( (z[1:n] - mean(z[1:n]))^3 ) / ( sum( (z[1:n] - mean(z[1:n]))^2 ) )^(3/2) } )
    tmp <- unlist( tmp )
    tail( cbind( M3, tmp ) )
  

    n <- 10
    sum( (z[1:n] - mean(z[1:n]))^3 ) / ( sum( (z[1:n] - mean(z[1:n]))^2 ) )^(3/2)
    skewFUN(z[1:n]) * (n / n^(3/2))
    M3[n-1, "avg_boot"]  
    
    
    
    # Expected lp_parts^2
    e2n
    s <- unlist( lapply( n_vec, FUN = function(n) { var(z[1:n]) } ) )
    e2n - 1/n_vec
    cbind( SE[ ,c("avg_boot", "exact_boot")]^2, s * ( e2n - 1/n_vec ) )
    
    
    
    
  }
  
  
  
  
 
  
  
  
  
  
  
  
  
  