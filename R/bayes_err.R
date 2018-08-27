bayeserr <- function (da, pi0, mu, Sigma, idx) {
  # bayes_error of using true parameters
  
  N <- nrow(da); 
  K <- length(pi0)
  idxbayes <- rep(NA, N)
  tempP <- matrix(0, N, K)
  if (length(dim(Sigma))==2) {
      for (j in 1:K){
        tempP[, j] <- pi0[j]*dmvnorm(da, mean = mu[, j], sigma = Sigma)
      }
  }else {
    for(j in 1:K){
      tempP[, j] <- pi0[j]*dmvnorm(da, mean=mu[, j], sigma=Sigma[, , j])
    }
  }
  idxbayes <- apply(tempP, 1, which.max)
  sum(idxbayes != idx)/N
}