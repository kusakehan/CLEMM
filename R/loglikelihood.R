Mixturell <- function(K, da, pi0, mu, Sigma) {
  #calculate observed loglikelihood
  
  N <- nrow(da)
  tempP <- matrix(0, N, K)
  if(length(dim(Sigma))==2) {
      for(ii in 1:K){
        tempP[, ii] <- dmvnorm(da, mean=mu[, ii], sigma=Sigma)
      }
  }else {
      for(ii in 1:K){
       tempP[, ii] <- dmvnorm(da, mean=mu[, ii], sigma=Sigma[, , ii])
      }
  }
  l <- sum(log(tempP %*% pi0))
  return(l)
}