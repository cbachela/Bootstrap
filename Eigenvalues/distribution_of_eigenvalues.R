  
  
  ############################################################################
  ### DISTRIBUTION OF SAMPLE EIGENVALUES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     05.07.2022
  # First version:    05.07.2022
  # --------------------------------------------------------------------------
  
  
  require(slolz)
  require(GPO)

  
  X <- getMSCIData( universe = "dm" )
  X <- X[ ,1:5]
  
  covmat <- cov(X)
  # covmat <- cor(X)
  E <- eigen(covmat)
  evalues <- E$values
  PCA <- pca( X, pcaCtrl(method = "eigen") )
  S <- getter( PCA, "sources" )
  Delta <- t(E$vectors) %*% covmat %*% E$vectors
  
  E$vectors %*% diag(E$values) %*% t(E$vectors) - covmat
  cov( X %*% E$vectors )
  cbind(diag( cov( X %*% E$vectors ) ), E$values)
  headleft(Delta)
  headleft(E$vectors)
  headleft(getter(PCA, "torsion"))
  torsion(covmat = covmat, method = "pc")
  
  plot( density( evalues ) )
  barplot( evalues )
  
  
  
  # Bootstrap in PC space
  
  n_sim <- 10^3
  n <- nrow(X)
  alpha <- rep(1, n)
  samples <- rdirichlet( n = n_sim, alpha = alpha )
  
  evalues_boot <- apply( samples, 1, function(w) { t(S^2) %*% w } )  
  plot( density(evalues_boot) )
  
  evalues_boot_avg <- apply( evalues_boot, 1, mean )  
    
  plot( density( evalues ) )
  lines( density( evalues_boot_avg ), col = 2 )
  lines( density( evalues_boot), col = 3 )
  
  
  boxplot( as.data.frame(t(evalues_boot)) )
  points( evalues, pch = 19, col = 2 )
  points( evalues_boot_avg, pch = 19, col = 3 )
  points( apply(tmp$boot.eigs, 2, mean), pch = 19, col = 4 )
  
  
  
  
  # remotes::install_github("HerveAbdi/data4PCCAR")
  require(data4PCCAR)
  boot.eigen
  
  #// re-computes an eigen transformation at each iteration. 
  #// I.e., eigenvectors vary.
  tmp <- boot.eigen( X = X, nIter = n_sim, center = TRUE )
  
  tmp$fixed.eigs
  
  
  plot( density( evalues ) )
  lines( density( evalues_boot_avg ), col = 2 )
  lines( density( evalues_boot ), col = 3 ) 
  lines( density( tmp$boot.eigs.sorted ), col = 4 )
  
  
  
  
  
  
  
  
  
