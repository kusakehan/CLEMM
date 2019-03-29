clemm_em <- function (da, K, u, iter, stopping=1e-7, opts=NULL, init, typ="G", add_eye="False") {
  
  # da is the data for clustering
  #  K is the number of clusters selected
  #  u is the number of envlope dimension selected
  #  iter is the maximum iteration
  #  stopping is the covergence criterion 
  #  opts is for 1D algorithm
  #  init is the list initial value for EM:
  #      init$centers; init$wt; init$cov
  #  typ indicates we use clemm or clemm-shared
  #      typ="G": general clemm
  #      typ="s": clemm shared with same covariance matrix
  #               across all clusters
  #  add_eye indicates whether add 0.01*diag(r) in clemm: True or False
  
  N <- dim(da)[1]; r <- dim(da)[2]
  muEst <- muEM <- array(NA, c(K, r, iter+1))
  a <- matrix(NA, K, iter+1)
  muEst[, , 1] <- t(init$centers)
  a[, 1] <- init$wt
  Sx <- (N-1)*cov(da)/N
  Xbar <- apply(da, 2, mean)
  ll <- rep(NA, iter)
  
  if (typ == "G") {
    if (is.null(opts$xtol))
      opts$xtol = 1e-10 else if (opts$xtol < 0 || opts$xtol > 1)
      opts$xtol = 1e-10 
      
      
    if (is.null(opts$gtol))
      opts$gtol = 1e-10 else if (opts$gtol < 0 || opts$gtol > 1)
      opts$gtol = 1e-10 
    
    if (is.null(opts$ftol))
      opts$ftol = 1e-10 else if (opts$ftol < 0 || opts$ftol > 1)
      opts$ftol = 1e-10 
          
    if (is.null(opts$mxitr))
      opts$mxitr = 2000
    
    if (is.null(opts$record))
      opts$record = 0
    
    SigmaEst <- array(NA, c(r, r, K, iter+1))
    S_tmp <- array(NA, c(r, r, K))
    SigmaEst[, , , 1] <- init$cov
    p1 <- rep(NA, K); p <- matrix(NA, K, N)
    S <- array(0, c(r, r, K, iter+1))
    ll[1] <- Mixturell(K, da, a[, 1], t(muEst[, , 1]), SigmaEst[, , , 1]) 
    
    for(n in 1:iter){
      tmp <- matrix(0, r, K)
      tt1 <- tt2 <- matrix(0, N, K)
      
      for (j in 1:K) {   
        tt1[, j] <- a[j, n]*dmvnorm(da, mean = muEst[j, , n], sigma = SigmaEst[, , j, n])
        for (s in 1:K) {
          tt2[, s] <- a[s, n]*dmvnorm(da, mean = muEst[s, , n], sigma = SigmaEst[, , s, n])
        }
        P <- apply(tt2, 1, sum)
        
        for(i in 1:N){
          p[j, i] <- tt1[i, j]/P[i]
        }
        
        tmp1 <- apply(p, 1, sum)[j]
        
        
        for (i in 1:N) {
          tmp[, j] <- p[j, i]*da[i, ] + tmp[, j]
        }
        
        
        muEM[j, , n] <- tmp[, j]/tmp1
        
        tmp3 <- matrix(0, r, r)
        for (i in 1:N) {
          tmp3 <- p[j, i]*((da[i, ]-muEM[j, , n]) %*% t((da[i, ]-muEM[j, , n])))+tmp3
        }
        
        a[j, n+1] <- tmp1/N
        S[, , j, n] <- tmp3/tmp1
        if(add_eye == "True"){S_tmp[, , j] <- S[, , j, n] + 0.01*diag(r)}else {S_tmp[, , j] <- S[, , j, n]}
      }
      
      Gammaest_init <- OptimballGBB1D(p, Sx, S_tmp, u, opts=NULL)
      Gammaest <- OptStiefelGBB(Gammaest_init, opts, FGfun, p, Sx, S_tmp)$X
      
      
      Gammaest0 <- qr.Q(qr(Gammaest), complete = TRUE)[, (u+1):r]
      
      
      
      for (j in 1:K){
        muEst[j, , n+1] <- Xbar + Gammaest %*% t(Gammaest) %*% matrix(muEM[j, , n]-Xbar)
        
        SigmaEst[, , j, n+1] <- Gammaest %*% (t(Gammaest) %*% S[, , j, n] %*% Gammaest) %*% t(Gammaest) + 
          Gammaest0 %*% (t(Gammaest0) %*% Sx %*% Gammaest0) %*% t(Gammaest0)
      }
      
      
      ll[n+1] <- Mixturell(K, da, a[, n+1], t(muEst[, , n+1]), SigmaEst[, , , n+1])
      #print(n)
      if((n>1) && (ll[n+1]-ll[n]<stopping)) break;
    }
    if ((ll[n+1]-ll[n]) < 0){
      lst = list(ll_seq=ll, a=a[, n], muEst=muEst[, , n], SigmaEst=SigmaEst[, , , n], iteration=n, Gammaest=Gammaest)
    }else {
      lst = list(ll_seq=ll, a=a[, n+1], muEst=muEst[, , n+1], SigmaEst=SigmaEst[, , , n+1], iteration=n+1, Gammaest=Gammaest)
    }
  }else if(typ == "S") {
    SigmaEst <- array(NA, c(r, r, iter+1))
    SigmaEst[, , 1] <- init$cov
    p1 <- rep(NA, K); p <- matrix(NA, K, N)
    S <- array(0, c(r, r, iter+1))
    ll[1] <- Mixturell(K, da, a[, 1], t(muEst[, , 1]), SigmaEst[, , 1]) 
    
    for(n in 1:iter){
      tmp <- matrix(0, r, K)
      tmp3 <- matrix(0, r, r)
      
      tt1 <- tt2 <- matrix(0, N, K)
      
      for (j in 1:K) {   
        tt1[, j] <- a[j, n]*dmvnorm(da, mean = muEst[j, , n], sigma = SigmaEst[, , n])
        for (s in 1:K) {
          tt2[, s] <- a[s, n]*dmvnorm(da, mean = muEst[s, , n], sigma = SigmaEst[, ,  n])
        }
        P <- apply(tt2, 1, sum)
        
        for(i in 1:N){
          p[j, i] <- tt1[i, j]/P[i]
        }
        
        tmp1 <- 0
        
        for (i in 1:N) {
          tmp1 <- p[j, i] + tmp1
          tmp[, j] <- p[j, i]*da[i, ]+tmp[, j]
        }
        
        
        muEM[j, , n] <- tmp[, j]/tmp1
        a[j, n+1] <- tmp1/N
        
        for (i in 1:N) {
          tmp3 <- p[j, i]*((da[i, ]-muEM[j, , n]) %*% t((da[i, ]-muEM[j, , n])))+tmp3
        }
      }
      
      S[, , n] <- tmp3/N
      Gammaest <- ECD(S[, , n], Sx, u)
      Gammaest0 <- qr.Q(qr(Gammaest), complete = TRUE)[, (u+1):r]
      
      
      
      for (j in 1:K){
        muEst[j, , n+1] <- Xbar + Gammaest %*% t(Gammaest) %*% matrix(muEM[j, , n]-Xbar)
      }
      
      SigmaEst[, , n+1] <- Gammaest %*% (t(Gammaest) %*% S[, , n] %*% Gammaest) %*% t(Gammaest) + 
        Gammaest0 %*% (t(Gammaest0) %*% Sx %*% Gammaest0) %*% t(Gammaest0)
      
      
      ll[n+1] <- Mixturell(K, da, a[, n+1], t(muEst[, , n+1]), SigmaEst[, , n+1])
      
      #print(n)
      if((n>1) && (ll[n+1]-ll[n]<stopping)) break;
    }
    if ((ll[n+1]-ll[n]) < 0){
      lst = list(ll_seq=ll, a=a[, n], muEst=muEst[, , n], SigmaEst=SigmaEst[, , n], iteration=n, Gammaest=Gammaest)
    }else {
      lst = list(ll_seq=ll, a=a[, n+1], muEst=muEst[, , n+1], SigmaEst=SigmaEst[, , n+1], iteration=n+1, Gammaest=Gammaest)
    }
  }
  return(lst)
}
