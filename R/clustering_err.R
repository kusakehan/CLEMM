clustering_err <- function(K, da, em_res, pi0=NULL, mu=NULL, Sigma=NULL, idx) {
  # da is data for clustering
  # pi0, mu, Sigma are true parameters
  # idx is the true label
  # em_res is a list contains the estimation sequence of em algorithm
  muEst = em_res$muEst
  SigmaEst = em_res$SigmaEst
  a = em_res$a
  
  tp <- permn(K);
  tp2 <- rep(NA, length(tp))
  N <- nrow(da);
  idxpred <- matrix(NA, N, length(tp))
  muer <- sigmaer <- pier <- 0
 
  tempP <- matrix(0, N, K)
  for(j in 1:K){
    if (length(dim(SigmaEst))==3) {
      tempP[, j] <- a[j]*dmvnorm(da, mean=muEst[j, ], sigma=SigmaEst[, , j])
    }else{tempP[, j] <- a[j]*dmvnorm(da, mean=muEst[j, ], sigma=SigmaEst)}
  }
  tmp <- apply(tempP, 1, which.max)
  for (t2 in 1:length(tp)) {
    
    idxpred[, t2] <- (tp[[t2]])[tmp]
    tp2[t2] <- sum(idxpred[, t2] != idx)/N
  }
  prederr <- min(tp2)
  
  ix <- which.min(tp2);
  od <- tp[[ix]]
  if (length(dim(Sigma))==3) {
    for (j in 1:K){
      muer <- norm(matrix(muEst[od[j], ]-mu[, j]), type = "F") + muer
      pier <- abs(a[od[j]]-pi0[j])+pier
      sigmaer <- norm((SigmaEst[, , od[j]]-Sigma[, , j]), type = "F")+sigmaer
    }
  }else if(!is.null(Sigma)){
    sigmaer <- norm((SigmaEst-Sigma), type = "F")
    for (j in 1:K){
      muer <- norm(matrix(muEst[od[j], ]-mu[, j]), type = "F") + muer
      pier <- abs(a[od[j]]-pi0[j])+pier
    }
  }
  return(list(cluster_err=prederr, mean_err=muer, wt_err=pier, cov_err=sigmaer, idxenv=idxpred[, ix]))
}
