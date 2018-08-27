#source("OptManiMulitBallGBB.R")

##################################################
#    get initial value for 1D algorithm          #
##################################################
get_ini1D <- function(p1, Sx, Sk) {
  r <- dim(Sx)[2]
  K <- dim(Sk)[3]
  v <- eigen(Sx)$vectors
  P_k <- apply(p1, 1, sum); n<- dim(p1)[2]
  
  for (k in 1:K) {
    v2 <- eigen(Sk[ , , k])$vectors
    v <- cbind(v2, v)
  }
  
  
  W0 <- v[, 1]
  Fw0 <- n*log(t(W0) %*% solve(Sx) %*% W0) 
  
  for (j in 1:K) {
    Fw0 <- Fw0 + P_k[j]*log((t(W0) %*% Sk[, , j] %*% W0)) 
  }
  
  
  for (i in 2:dim(v)[2]) {
    W <- v[, i]
    Fw <- n*log(det(t(W) %*% solve(Sx) %*% W))
    for (j in 1:K) {
      Fw <- Fw + P_k[j]*log((t(W) %*% Sk[, , j] %*% W))
    }
    if (Fw < Fw0){
      W0 <- W
      Fw0 <- Fw
    }
  }
  return(W0)
}
##################################################
#         1D optimization function               #
##################################################
fun1D <- function(W, p1, Sx, Sk) {
  r <- dim(Sx)[2]
  K <- dim(Sk)[3]
  P_k <- apply(p1, 1, sum); n<- dim(p1)[2]
  
  f <- n*log(t(W) %*% solve(Sx) %*% W) 
  for (k in 1:K){
     f <- f+P_k[k]*log(t(W) %*% Sk[, , k] %*% W) 
  }
  df <- n*solve(Sx) %*% W %*% solve(t(W) %*% solve(Sx) %*% W) 
  for (k in 1:K){
    df <- df + P_k[k]*Sk[, , k] %*% W %*% solve(t(W) %*% Sk[, , k] %*% W)
  }
  df <- 2*df
  return(list(F = f, G = df))
}


##################################################
#         1D optimization solve for gamma        #
##################################################
ballGBB1D <- function(p1, Sx, Sk, opts=NULL) {
  W0 <- get_ini1D(p1, Sx, Sk)
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
        
  X <- OptManiMulitBallGBB(W0, opts, fun1D, p1, Sx, Sk)$X
  return(X)
}

##################################################
#  1D optimization solve for envelope basis      #
##################################################   

OptimballGBB1D <- function(p1, Sx, Sk, d, opts=NULL) {
  
  if(dim(Sx)[1]!=dim(Sx)[2]){
    {Sx = Sx %*% t(Sx)}
  }
  r = dim(Sx)[2]
  if(d < r){
    Sknew <- Sk
    Sxnew <- Sx
    G <- matrix(0, r, d)
    G0 <- diag(r)
    for(k in 1:d){
      gk <- ballGBB1D(p1, Sxnew, Sknew, opts)
      G[, k] <- G0 %*% gk
      G0 <- qr.Q(qr(G[, 1:k]),complete=T)[,(k+1):r]
      Sxnew <- t(G0) %*% Sx %*% G0
      Sknew <- array(0, c((dim(Sknew)[1]-1), (dim(Sknew)[2]-1), K))
      for (j in 1:K) {
        Sknew[ , ,j] <- t(G0) %*% Sk[ , , j] %*% G0
      }
    }
    Ghat <- G
  } else { Ghat <- diag(r) }
  return(Ghat)
}