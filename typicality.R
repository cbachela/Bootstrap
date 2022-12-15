  
  
  ############################################################################
  ### BOOTSTRAP - TYPICALITY
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     24.04.2021
  # First version:    16.04.2021
  # --------------------------------------------------------------------------
  
  
  # Toy example, k = 2
  # Toy example, k = 3
  # Typical partitions, n = 6 (partitions form the alphabet)
  # Typical partitions, n = 80
  # n = 6
  # Loop over n and B, check if AEP holds
  
  
  
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
  
  
  # --------------------------------------------------------------------------
  # FUNCTIONS
  # --------------------------------------------------------------------------
  
  # Sequency occurrence probability
  probFUN <- function( x, prob )
  {
    ans <- 1
    for ( i in seq(along = prob) ) {
      ans <- ans * prob[i]^(sum(x == i))
    }
    return( ans )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Toy example, k = 2
  # --------------------------------------------------------------------------
  
  # See https://en.wikipedia.org/wiki/Typical_set
  
  # See https://mc-stan.org/users/documentation/case-studies/curse-dims.html
  
  # Suppose we have a binary sequence of N elements, z = z_i,.,z_n, with zn???{0,1} and let
  # y = \sum_{i=1}^n = z_i.
  
  # The repeated Bernoulli trial probability of z with a chance of success ?????[0,1] is 
  # given by the probability mass function
  
  # p(z?????) = \prod_{i=1}^n Bernoulli(z_i?????) = ??^y (1?????)^{(N???y)}.
  
  # The number of ways in which a binary sequence with a total of y successes 
  # out of n trials can arise is the total number of ways a subset of y elements 
  # may be chosen out of a set of n elements, which is given by the binomial coefficient,
  
  # $\choose(n, y) = \frac{n!}{y! (n???y)!}$.
  
  
  
  set.seed(1234)
  n <- 10^2
  B <- 10^2 + 1
  prob <- 0.1
  
  # Entropy of source
  H <- entropy( p = c(prob, 1-prob), 
                exponential = FALSE, 
                base = 2 )
  H
  
  # Most likely sequence
  y <- rep(1, n)
  p_y <- (1 - prob)^n
  p_y
  
  # Most likely sequence does not belong to the typical set
  # because
  -1/n * log(p_y, base = 2)
  H
  
  # and because its average logarithmic probability cannot come arbitrarily 
  # close to the entropy of the source no matter how large we take n.
  a <- -1/n * log( (1 - prob)^(1:10^4), base = 2 )
  plot( a )
  abline( h = H )
  
  # Sample
  samples <- matrix( NA, nrow = B, ncol = n )
  for ( i in 1:B ) {
    x_unif <- runif( n )
    samples[i, ] <- as.numeric( x_unif > prob )
  }
  
  # Sequence occurance probabilities (Bernoulli probabilities)
  p_k <- apply( samples, 1, function(x) { (1 - prob)^sum(x) * prob^(n - sum(x)) } )
  p_k  
  
  # Check typicality
  plot( -1/n * log(p_k, base = 2) )
  abline( h = H, col = 2 )  
  
  # Or
  2^(-n * H)
  plot( p_k )
  abline( h = 2^(-n * H), col = 2 )
  
  quantile( -1/n * log(p_k, base = 2), c(0.1, 0.9) )
  H
  

  
  # Link to binomial distribution
  
  k <- 0:n
  b <- choose(n, k) # Number of members per leave / permutation per partition (binomial coefficient)
  d <- (1 - prob)^k * prob^(n - k) # Bernoulli probability
  
  plot( x = k, y = b, type = "o" )
  plot( x = k, y = d, type = "o" )
  plot( x = k, y = b *  d, type = "o" )
  
  sum( d )
  sum( b * d )
  
  
  # Using dbinom
  
  dbinom( x = 0:n, size = n, prob = prob )
  sum( dbinom( x = 0:n, size = n, prob = prob ) )  # should be 1
  
  tmp <- apply( samples, 1, function(x) { dbinom( x = sum(x), size = n, prob = (1 - prob) ) } )
 
  head(tmp)
  head( p_k * apply( samples, 1, function(x) { choose(n, sum(x)) } ) )  # same same
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Toy example, k = 3
  # --------------------------------------------------------------------------
 
  
  # $X = (X_1, ..., X_k)$, $X_j$ the number of times event $E_k$ occurs,
  # is said to have a multinomial distribution with index 
  # $n$ and parameter $\theta$, i.e., $X \sim Mult(n, \theta)$.
  # The individual components of a multinomial random vector are binomial and 
  # have a binomial distribution,
  # $X_1 \sim Binom(n, \theta_1)$.
  
  
  
  set.seed(1234)
  n <- 3
  B <- 10^2 + 1
  k <- 3
  z <- 1:3
  prob <- c(0.1, 0.5, 0.4)
 
  # Partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat )
  prob <- draws2prob( draws = draws )
  prob_tot <- prob * draws[2, ]
  prob_dmult <- apply( parts_mat, 2, dmultinom, prob = rep(1/n, n), log = FALSE )
  
  # Entropy of source
  H <- entropy( p = prob,
                exponential = FALSE, 
                base = 2 )
  H
  
  # Sample
  samples <- matrix( sample( x = 1:k, size = B * n, prob = prob, replace = TRUE ),
                     nrow = B, ncol = n )
  
  # Checks
  sum( samples == 1 ) / length(samples)
  sum( samples == 2 ) / length(samples)
  sum( samples == 3 ) / length(samples)
  
  # Sequence occurrence probabilities
  p_k <- apply( samples, 1, probFUN, prob = prob )
  
  
  # Check typicality
  plot( -1/B * log(p_k, base = 2) )
  abline( h = H, col = 2 )  
  
  # Or
  plot( p_k )
  abline( h = 2^(-n * H), col = 2 )
  
  quantile( -1/B * log(p_k, base = 2), c(0.1, 0.9) )
  H
  
  
  
  
  tmp <- as.numeric( apply( samples, 1, dbinom, size = n, prob = prob ) )
  head(tmp)  
  head(p_k)
  
  
  
  
  # Link to multinomial distribution
  
  counts <- matrix( 0, nrow = B, ncol = k ) 
  for ( i in 1:nrow(counts) ) {
    tbl <- table(samples[i, ])
    counts[i, as.numeric(names(tbl))] <- tbl
  }
  head(counts)
  p_dmult <- apply( counts, 1, dmultinom, size = n, prob = prob )
  p_dmult
  
  
  mcoeff111 <- factorial(n) / prod( factorial(c(1,1,1)) )
  mcoeff210 <- factorial(n) / prod( factorial(c(1,2,0)) )
  mcoeff300 <- factorial(n) / prod( factorial(c(3,0,0)) )
  mcoeff111
  mcoeff210
  mcoeff300
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  dmultinom( c(0,2,1), size = n, prob = rep(1/k, k) )
  dmultinom( c(2,0,1), size = n, prob = rep(1/k, k) ) * 6
  dmultinom( c(1,1,1), size = n, prob = rep(1/k, k) )
  dmultinom( c(3,0,0), size = n, prob = rep(1/k, k) ) * 3
  
  ###
  # prob <- rep(1/k, k)
  probFUN( x = c(1,1,3), prob = prob )
  probFUN( x = c(1,1,3), prob = prob ) * factorial(n) / prod( factorial(c(2,0,1)) )
  dmultinom( c(2,0,1), size = n, prob = prob )
  
  probFUN( x = c(1,1,1), prob = prob )
  probFUN( x = c(1,1,1), prob = prob ) * factorial(n) / prod( factorial(c(3,0,0)) )
  dmultinom( c(3,0,0), size = n, prob = prob )
  dmultinom( c(0,3,0), size = n, prob = prob )
  probFUN( x = c(2,2,2), prob = prob ) * factorial(n) / prod( factorial(c(0,3,0)) )
  
  
  probFUN( x = c(1,2,3), prob = prob )
  probFUN( x = c(1,2,3), prob = prob ) * factorial(n) / prod( factorial(c(1,1,1)) )
  dmultinom( c(1,1,1), size = n, prob = prob )
  ###
  
  
  
  
  
  
  x <- as.numeric(counts[1, ])
  coeff <- factorial(n) / prod( factorial(x) )
  coeff * p_k
  
  
  # Or, more clearly
  plot( p_dmult )
  abline( h = 2^(-n * H), col = 2 )
  
  head( cbind( p_k, p_dmult ) ) 
    
    
  
  
  # --------------------------------------------------------------------------
  # Typical partitions, n = 6 (partitions form the alphabet)
  # --------------------------------------------------------------------------
  
  n <- 6
  B <- 80
  M <- 10^2
  
  # Enumerate all partitions
  parts_mat <- partitions::parts( n )
  draws <- gridDraws( parts_mat = parts_mat )
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
  

  # Entropy of probability distribution of all partitions
  H <- entropy( p = prob_tot, exponential = FALSE, base = 2, eps = 1e-15 )
  H  
  sum( prob_tot * log( 1 / prob_tot, base = 2 ) )

  
  # Run bootstrap procedure m times
  # and compute probability to observe the sequence of B partitions
  p_m <- numeric(M)
  for ( m in 1:M ) {
    
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample(1:n, replace = TRUE)
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl
    }
    gridnodes_sort <- t( apply( gridnodes, 1, function(x) { rev(sort(x)) } ) )  
    
    # Match sampled gridnodes to analytic level sets, aka partitions
    id_match <- apply( gridnodes_sort, 1, prodlim::row.match, 
                       table = t(as.matrix(parts_mat[ ,1:min(ncol(parts_mat), 10^4)])) )
    
    p_m[m] <- probFUN( x = id_match, prob = prob_tot )
  }
  
  
  # Check typicality
  plot( -1/B * log(p_m, base = 2) )
  abline( h = H, col = 2 )  
  
  # Or
  plot( p_m )
  abline( h = 2^(-B * H), col = 2 )
  
  # quantile( -1/B * log(p_m, base = 2), c(0.01, 0.99) )
  range( -1/B * log(p_m, 2) )
  H
  
  range( log(p_m, base = 2) )
  2^(B * H)
  
  
  
  # Most likely sequence does not belong to the typical set
  p_ml <- max(prob_tot^B)
  log(p_ml, 2)
  range( log(p_m, 2))
  H
  
  -1/B * log(p_ml, 2)
  H
  
  

    
  # --------------------------------------------------------------------------
  # Typical partitions, n = 80
  # --------------------------------------------------------------------------
  
  n <- 80
  B <- 10^2
  M <- 10^2
  
  env <- readRDS( file = paste0(wd, "waRehouse/bootstrap_grid_lp_norm_n80.rds") )
  n <- env$n
  B <- env$B
  gridnodes <- env$gridnodes
  parts_mat <- env$parts_mat
  draws <- env$draws
  prob <- env$prob
  prob_dmult <- env$prob_dmult
  prob_tot <- env$prob_tot
  
  
  
  # Entropy of probability distribution of all partitions
  H <- entropy( p = prob_tot, exponential = FALSE, base = 2, eps = 1e-15 )
  H  
  sum( prob_tot * log( 1 / prob_tot, base = 2 ) )
  
  
  # Sample a sequence of partitions
  # 1) via ordinary bootstrap
  
  gridnodes <- matrix(0, nrow = B, ncol = n)
  for ( i in 1:B ) {
    idx <- sample(1:n, replace = TRUE)
    tbl <- table(idx)
    gridnodes[i, as.numeric(names(tbl))] <- tbl
  }
  gridnodes_sort <- t( apply( gridnodes, 1, function(x) { rev(sort(x)) } ) )  
  
  # Match sampled gridnodes to analytic level sets, aka partitions
  id_match <- apply( gridnodes_sort, 1, prodlim::row.match, 
                     table = t(as.matrix(parts_mat[ ,1:min(ncol(parts_mat), 10^4)])) )
  
  # 2) via direct sampling from parts_mat with prob = prob_tot
  idx <- sample( 1:ncol(parts_mat), size = B, replace = TRUE, prob = prob_tot )
  p_idx <- probFUN( x = idx, prob = prob_tot )
  p_id <- probFUN( x = id_match, prob = prob_tot )
  
  
  
  # Check typicality
  plot( -1/B * log(p_id, base = 2) )
  abline( h = H, col = 2 )  
  
  # Or
  plot( p_k )
  abline( h = 2^(-n * H), col = 2 )
  
  quantile( -1/B * log(p_id, base = 2), c(0.01, 0.99) )
  H
  
  
  
  ###
  
  p_m <- numeric(M)
  
  for ( m in 1:M ) {
    
    gridnodes <- matrix(0, nrow = B, ncol = n)
    for ( i in 1:B ) {
      idx <- sample(1:n, replace = TRUE)
      tbl <- table(idx)
      gridnodes[i, as.numeric(names(tbl))] <- tbl
    }
    gridnodes_sort <- t( apply( gridnodes, 1, function(x) { rev(sort(x)) } ) )  
    
    # Match sampled gridnodes to analytic level sets, aka partitions
    id_match <- apply( gridnodes_sort, 1, prodlim::row.match, 
                       table = t(as.matrix(parts_mat[ ,1:min(ncol(parts_mat), 10^4)])) )
    
    p_m[m] <- probFUN( x = id_match, prob = prob_tot )
    
  }
  
  
  # Check typicality
  plot( -1/B * log(p_m, base = 2) )
  abline( h = H, col = 2 )  
  
  # Or
  plot( p_m )
  abline( h = 2^(-B * H), col = 2 )
  
  quantile( -1/B * log(p_m, base = 2), c(0.01, 0.99) )
  H
  
  range(log(p_m, base = 2))
  2^(B * H)
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # n = 6
  # --------------------------------------------------------------------------
  
  n <- 6
  prob <- rep(1/n, n)
  x <- c(1, 1, 2, 2, 3, 4)
  counts <- c(2, 2, 1, 1, 0, 0)
  table(x)
  p <- (1/n)^n
  probFUN( x = x, prob = prob )
  p
  num <- factorial( sum(counts) )
  tbl <- table(counts)
  a <- num / prod( factorial( counts ) )    # multinomial coefficient!
  b <- num / prod( factorial( tbl ) )
  p * b
  p * a
  dmultinom( counts, size = n, prob = prob )
  
  
  ddirichlet( x = x / n, alpha = rep(1, n) )
  prod( (x / n)^prob )

  
 
  
  
  # --------------------------------------------------------------------------
  # Loop over n and B, check if AEP holds
  # --------------------------------------------------------------------------
  
  # Data
  n_vec <- 3:40
  b_vec <- seq(from = 10, to = 10^3, length.out = 100)
  z <- rt( df = 4, n = max(n_vec) )
  
  for ( i in seq(along = n_vec) ) {
    for ( j in seq(along = b_vec ) ) {
      
      n <- n_vec[i]
      b <- b_vec[j]
      
      # Partitions
      parts_mat <- partitions::parts( n )
      draws <- gridDraws( parts_mat = parts_mat, cnames = FALSE, mpfr = FALSE )
      prob <- draws2prob( draws = draws )
      prob_tot <- prob * draws[2, ]
      
      # Entropy
      H <- entropy( p = prob_tot, exponential = FALSE, base = 2, eps = 1e-15 )
      
      # Bootstrap
      gridnodes <- matrix(0, nrow = B, ncol = n)
      for ( i in 1:B ) {
        idx <- sample(1:n, replace = TRUE)
        tbl <- table(idx)
        gridnodes[i, as.numeric(names(tbl))] <- tbl
      }
      gridnodes_sort <- t( apply( gridnodes, 1, function(x) { rev(sort(x)) } ) )  
      
      # Match sampled gridnodes to analytic level sets, aka partitions
      id_match <- apply( gridnodes_sort, 1, prodlim::row.match, 
                         table = t(as.matrix(parts_mat[ ,1:min(ncol(parts_mat), 10^4)])) )
      
      p_m <- probFUN( x = id_match, prob = prob_tot )
      H
      -1/B * log(p_m, 2)
      
      
    }
  }
  
  
  
  ################################################### 
  
  # Plot 
  require(rgl)
  
  n <- 4
  n_sim <- 10^7
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  eps <- 1e-03
  idx <- apply( samples, 1, function(x) { abs( sum(x^2) - 2/(n+1) ) < eps } )
  
  
  # m <- 10^5 * 2
  # samples_fev <- fevBias( samples[1:m, ], q = 100 )
  # rgl::plot3d( x = samples[1:m ,1], y = samples[1:m ,2], z = samples[1:m ,3], col = "white",
  #              size = 1, xlab = "", ylab = "", zlab = "", axes = FALSE )
  # points3d(  x = samples[idx ,1], y = samples[idx ,2], z = samples[idx ,3], col = "orange" )
  # points3d(  x = samples_fev[ ,1], y = samples_fev[ ,2], z = samples_fev[ ,3], col = "grey" )
  
  
  
  
  m <- 10^5 * 2
  samples_fev <- fevBias( samples[1:m, ], q = 100 )
  rgl::plot3d( x = samples_fev[ ,1], y = samples_fev[ ,2], z = samples_fev[ ,3], col = "grey",
               size = 2, xlab = "", ylab = "", zlab = "", axes = FALSE )
  points3d(  x = samples[idx ,1], y = samples[idx ,2], z = samples[idx ,3], size = 0.7, col = "orange" )

  
  
  
  
  
  
    