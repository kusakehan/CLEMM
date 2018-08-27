##################################################
#         ECD objective function                 #
##################################################
objF_ECD <- function(A, B, w) {
  fk <- log(t(w) %*% A %*% w) + log(t(w) %*% B%*% w) - 2*log(t(w) %*% w)
  return(fk)
}

##################################################
#    get initial value for ECD algorithm         #
##################################################
ECDini <- function(S, Sx) {
  p <- dim(S)[2]
  eigS <- eigen(S); eigSx <- eigen(Sx)
  v1 <- eigS$vectors 
  v2 <- eigSx$vectors 
  v <- cbind(v1, v2)
  W0 <- v[, 1]
  Fw0 <- log(t(W0) %*% solve(Sx) %*% W0) + log(t(W0) %*% S %*% W0)
  for (i in 2:(2*p)) {
    W <- v[, i]
    Fw <- log(t(W) %*% solve(Sx) %*% W) + log(t(W) %*% S %*% W)
    if (Fw < Fw0){
      W0 <- W
      Fw0 <- Fw
    }
  }
  return(W0)
}

##################################################
#     ECD algorithm for solving fk               #
##################################################
optimeECD <- function(A, B, w0, maxiter) {
  p <- length(w0)
  epsilon <- 1e-08
  eigA <- eigen (A + t(A)); 
  ## already in descending order###
  Gp <- eigA$vectors; dn <- eigA$values 
  dn <- diag(dn/2)
  v0 <- t(Gp) %*% w0
  GBG <- t(Gp) %*% B %*% Gp
  fk <- objF_ECD(dn, GBG, v0)
  v <- v0
  for (iter in 1:maxiter) {
    flg <- 0
    alpha <- 1/(t(v) %*% dn %*% v)
    beta <- 1/(t(v) %*% GBG %*% v)
    delta <- 1/(t(v) %*% v)
    A1 <- as.numeric(alpha)*dn
    B1 <- as.numeric(beta)*GBG
    for (j in 1:p) {
      AB1 <- A1[j, j] + B1[j, j]
      if ((2*delta - AB1) != 0) {
        v[j] <- (t(v) %*% A1[, j] + t(v) %*% B1[, j] - AB1*v[j])/(2*delta - AB1)
        flg <- flg + 1
      }
      if (objF_ECD(dn, GBG, v) > (objF_ECD(dn, GBG, v0)) + epsilon) {
        v <- v0
        flg <- flg - 1
      }
    }
    fk1 <- objF_ECD(dn, GBG, v)
    if ((abs (fk - fk1)) < epsilon) break
    fk <- fk1
  }
  w <- Gp %*% v
  w <- w/as.numeric(sqrt(t(w) %*% w))
  return(w)
}

##################################################
#   estimating envelop subspace                  #
##################################################
ECD1st <- function (S, Sx, maxiter){
  gamma <- optimeECD(S, solve(Sx), ECDini(S, Sx), maxiter)
  return(gamma)
}

##################################################
#           ECD algorithm                        #
##################################################
ECD <- function(S, Sx, d, maxiter = 2000){
  # S>0, Sx>0 and is symmetric
  # dimension of the envelope is d
  # based on inv(Sx) and (S)
  p <- dim(S)[2]
  Snew <- S
  Sxnew <- Sx
  G <- matrix(0, p, d)
  G0 <- diag(p)
  for (k in 1:d) {
    gk <- ECD1st(Snew, Sxnew, maxiter)
    G[, k]<- G0 %*% gk
    G0 <- qr.Q(qr(G[, 1:k]),complete=T)[,(k+1):p]
    Snew <- t(G0) %*% S %*% G0
    Sxnew <- t(G0) %*% Sx %*% G0
  }
  Ghat <- G
  return(Ghat)
}

