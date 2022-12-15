


  n <- 10^2
  N <- 10^3
  set.seed(1234)
  X <- matrix( rnorm( n = n*N ), nrow = N, ncol = n )
  
  alpha_n <- rep(1, n)
  # alpha_n <- (1:n) / sum(1:n) * n
  w_bar <- alpha_n / n
  alpha_N <- rep(1, N)
  # alpha_N <- (1:N) / sum(1:N) * N
  w_bar_tilde <- alpha_N / N
  
  # Simulate portfolio returns
  sim_portf <- X %*% w_bar
  
  # Compute average asset returns
  mu <- apply( X, 2, mean )
  
  sum( mu * w_bar ); mean( sim_portf )

  
  # Sample weights from Simplex
  P <- rdirichlet( n = 10^4, alpha = alpha_n )
  apply(P, 1, sum)
  theta <- apply( P, 1, function(p) { sum(p * mu) } )
  
  P_tilde <- rdirichlet( n = 10^4, alpha = alpha_N )
  theta_2  <- apply( P_tilde, 1, function(p) { sum(p * sim_portf) } )

  
  plot( density( theta) )
  lines( density(theta_2), col = 2 )
  
  
  plot( density( X ) )
  lines( density( theta) )
  lines( density(theta_2), col = 2 )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
  
