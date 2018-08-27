data_generation <- function (r, u, N, K, md, seed=999) {
  
  # simulated data genration
  # r is dimension of data, u is the envelope dimension
  # N is number of data to be generated
  # K is number of clusters
  # md is type of model in the paper, choose from 1 to 5
  # seed is the random seed for generation process
  
  set.seed(seed)
  Gamma <- matrix(runif(r*u), r, u)
  Gamma <- qr.Q(qr(Gamma))
  Gamma0 <- qr.Q(qr(Gamma), complete = T)[, (u+1):r]

  
  if (md == 1) {
    Omega <- array(NA, c(u, u, K))
    for (j in 1:K) {
      tmp_omega <- matrix(runif(u^2), u, u)
      tmp_omega <- tmp_omega %*% t(tmp_omega)
      Omega[, , j] <- tmp_omega/norm(tmp_omega, type = "F")
    }
    
    Omega0 <- matrix(runif((r-u)^2), r-u, r-u)
    Omega0 <- Omega0 %*% t(Omega0)
    Omega0 <- Omega0/norm(Omega0, type = "F")
    
    Sigma <- array(NA, c(r, r, K))
    for (j in 1:K) {
      Sigma[, , j] <- exp(-j)*Gamma %*% Omega[, , j] %*% t(Gamma) + 
        Gamma0 %*% Omega0 %*% t(Gamma0)
      Sigma[, , j] <- 1.25*Sigma[, , j]/norm(Sigma[, , j], type = "F")
    }
  }else if(md == 2) {
      Omega <- array(NA, c(u, u, K))
      for (j in 1:K) {
        tmp_omega <- matrix(runif(u^2), u, u)
        tmp_omega <- tmp_omega %*% t(tmp_omega)
        Omega[, , j] <- tmp_omega/norm(tmp_omega, type = "F")
      }
      
      Omega0 <- matrix(runif((r-u)^2), r-u, r-u)
      Omega0 <- Omega0 %*% t(Omega0)
      Omega0 <- Omega0/norm(Omega0, type = "F")
      
      Sigma <- array(NA, c(r, r, K))
      for (j in 1:K) {
        Sigma[, , j] <- exp(-0.1*j)*Gamma %*% Omega[, , j] %*% t(Gamma) + 
          Gamma0 %*% Omega0 %*% t(Gamma0)
        Sigma[, , j] <- 3.5*Sigma[, , j]/norm(Sigma[, , j], type = "F")
      }
  }else if(md == 3) {
      Omega <- array(NA, c(u, u, K))
      for (j in 1:K) {
        tmp_omega <- matrix(runif(u^2), u, u)
        tmp_omega <- tmp_omega %*% t(tmp_omega)
        Omega[, , j] <- tmp_omega/norm(tmp_omega, type = "F")
      }
      
      Omega0 <- matrix(runif((r-u)^2), r-u, r-u)
      Omega0 <- Omega0 %*% t(Omega0)
      Omega0 <- Omega0/norm(Omega0, type = "F")
      
      Sigma <- array(NA, c(r, r, K))
      for (j in 1:K) {
        Sigma[, , j] <- Gamma %*% Omega[, , j] %*% t(Gamma) + 
          Gamma0 %*% Omega0 %*% t(Gamma0)
        Sigma[, , j] <- 12*Sigma[, , j]/norm(Sigma[, , j], type = "F")
      }
  }else if(md == 4){
      Omega <- matrix(runif(u^2), u, u)
      Omega <- Omega %*% t(Omega)
      Omega <- Omega/norm(Omega, type = "F")
      
      Omega0 <- matrix(runif((r-u)^2), r-u, r-u)
      Omega0 <- Omega0 %*% t(Omega0)
      Omega0 <- Omega0/norm(Omega0, type = "F")
    
      Sigma <- Gamma %*% Omega %*% t(Gamma) + 
        Gamma0 %*% Omega0 %*% t(Gamma0)
      Sigma <- 0.25*Sigma/norm(Sigma, type = "F")
  }else if(md == 5){
      Omega <- matrix(runif(u^2), u, u)
      Omega <- Omega %*% t(Omega)
      Omega <- Omega/norm(Omega, type = "F")
      
      Omega0 <- matrix(runif((r-u)^2), r-u, r-u)
      Omega0 <- Omega0 %*% t(Omega0)
      Omega0 <- Omega0/norm(Omega0, type = "F")
    
      Sigma <- Gamma %*% Omega %*% t(Gamma) + 
        Gamma0 %*% Omega0 %*% t(Gamma0)
      Sigma <- 2*Sigma/norm(Sigma, type = "F")
  }
  
  mu <- matrix(NA, r, K)
  eta <- matrix(NA, u, K)
  
  for (j in 1:K) {
    eta[, j] <- rnorm(u)
    mu[, j] <- Gamma %*% eta[, j]
  }
  
  return(list(mu=mu, Sigma=Sigma, Gamma=Gamma))
}