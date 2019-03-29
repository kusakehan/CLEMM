gmm_em <- function (da, K, iter, stopping=1e-7, init, typ="G", add_eye="False") {
  # da is the data for clustering
  #  K is the number of clusters selected
  #  iter is the maximum iteration
  #  stopping is the covergence criterion 
  #  init is the list initial value for EM:
  #      init$centers; init$wt; init$cov
  #  typ indicates we use gmm or gmm-shared
  #      typ="G": general gmm
  #      typ="S": gmm shared with same covariance matrix
  #               across all clusters
  #  add_eye indicates whether add 0.01*diag(r) in gmm: True or False
  
  N <- dim(da)[1]; r <- dim(da)[2]
  muEst <- muEM <- array(NA, c(K, r, iter+1))
  a <- matrix(NA, K, iter+1)
  muEst[, , 1] <- t(init$centers)
  a[, 1] <- init$wt
  Sx <- (N-1)*cov(da)/N
  Xbar <- apply(da, 2, mean)
  ll <- rep(NA, iter)
  
  if (typ == "G") {
    SigmaEst <- array(NA, c(r, r, K, iter+1))
    SigmaEst[, , , 1] <- init$cov
    p1 <- rep(NA, K); p <- matrix(NA, K, N)
    ll[1] <- Mixturell(K, da, a[, 1], t(muEst[, , 1]), SigmaEst[, , , 1]) 
    
    for(n in 1:iter){
      tmp <- matrix(0, r, K)
      tt1 <- tt2 <- matrix(0, N, K)
      
      for (j in 1:K) {   # 2 classes
        tt1[, j] <- a[j, n]*dmvnorm(da, mean = muEst[j, , n], sigma = SigmaEst[, , j, n])
        for (s in 1:K) {
          tt2[, s] <- a[s, n]*dmvnorm(da, mean = muEst[s, , n], sigma = SigmaEst[, , s, n])
        }
        P <- apply(tt2, 1, sum)
        
        for(i in 1:N){
          p[j, i] <- tt1[i, j]/P[i]
        }
        
        tmp1 <- 0
        
        
        
        tmp1 <- apply(p, 1, sum)[j]
        
        
        for (i in 1:N) {
          tmp[, j] <- p[j, i]*da[i, ] + tmp[, j]
        }
        
        
        muEst[j, , n+1] <- tmp[, j]/tmp1
        
        tmp3 <- matrix(0, r, r)
        for (i in 1:N) {
          tmp3 <- p[j, i]*((da[i, ]-muEst[j, , n+1]) %*% t((da[i, ]-muEst[j, , n+1])))+tmp3
        }
        
        a[j, n+1] <- tmp1/N
        SigmaEst[, , j, n+1] <- tmp3/tmp1
        if (add_eye == "True") {SigmaEst[, , j, n] <- SigmaEst[, , j, n] + 0.01*diag(r)}
      }
      ll[n+1] <- Mixturell(K, da, a[, n+1], t(muEst[, , n+1]), SigmaEst[, , , n+1])
      
      #print(c(n))
      if ((ll[n+1]-ll[n]) < stopping) break;
    }
    sig = SigmaEst[, , , n]
  }else if (typ == "S") {
    SigmaEst <- array(NA, c(r, r, iter+1))
    SigmaEst[, , 1] <- init$cov
    p1 <- rep(NA, K); p <- matrix(NA, K, N)
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
        
        
        muEst[j, , n+1] <- tmp[, j]/tmp1
        a[j, n+1] <- tmp1/N
        
        for (i in 1:N) {
          tmp3 <- p[j, i]*((da[i, ]-muEst[j, , n+1]) %*% t((da[i, ]-muEst[j, , n+1])))+tmp3
        }
      }
      
      SigmaEst[, , n+1] <- tmp3/N
      
      ll[n+1] <- Mixturell(K, da, a[, n+1], t(muEst[, , n+1]), SigmaEst[, , n+1])
      
      #print(c(n))
      if ((ll[n+1]-ll[n]) < stopping) break;
    }
    sig = SigmaEst[, , n]
  }
  return(list(a=a[, n], muEst=muEst[, , n], SigmaEst=sig, iteration=n, ll_seq=ll))
}
