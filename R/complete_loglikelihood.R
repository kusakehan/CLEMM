mixll_gamma <- function(p1, Gamma, M, U) {
  # calculate complete loglikelihood
  # approximated by Gamma function
  # p1 is mixing weights
  
  l <- log(det(t(Gamma) %*% solve(U) %*% Gamma))
  if (length(dim(M))==2) {
     l <- log(det(t(Gamma) %*% M %*% Gamma)) + l
  }else {
    K <- dim(M)[3]
    for (j in 1:K) {
      l <- l + p1[j]*log(det(t(Gamma) %*% M[, , j] %*% Gamma))
    }
  }
  return(l)
}

