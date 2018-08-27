FGfun <- function(W, p1, Sx, Sk) {
  r <- dim(Sx)[2]
  K <- dim(Sk)[3]
  P_k <- apply(p1, 1, sum); n<- dim(p1)[2]
  
  f <- n*log(det(t(W) %*% solve(Sx) %*% W) )
  for (k in 1:K){
    f <- f+P_k[k]*log(det(t(W) %*% Sk[, , k] %*% W))
  }
  df <- n*solve(Sx) %*% W %*% solve(t(W) %*% solve(Sx) %*% W) 
  for (k in 1:K){
    df <- df + P_k[k]*Sk[, , k] %*% W %*% solve(t(W) %*% Sk[, , k] %*% W)
  }
  df <- 2*df
  return(list(F = f, G = df))
}