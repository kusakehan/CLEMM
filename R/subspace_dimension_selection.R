env_dim_selection <- function(dim_rng, da, K, iter, stopping=1e-7, opts=NULL, init, typ="G", add_eye="False"){
  
  # dim_rng is the range of envelope dimension 
  # we want to choose from
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
  
  len <- length(dim_rng)
  N <- dim(da)[1]
  Sx <- (N-1)*cov(da)/N
  awe <- rep(NA, len)
  for (i in 1:len) {
   u0 <- dim_rng[i]
   res <- clemm_em(da, K, u0, iter, stopping, opts, init, typ, add_eye)
 
   Gammaest = res$Gammaest
   Sigmaest = res$SigmaEst
   a = res$a
   if (typ=="G"){
     df <- r+(K-1)*u0+r*(r+1)/2+(K-1)*u0*(u0+1)/2+(K-1)
   }else { df <- r +(K-1)*u0+r*(r+1)/2+(K-1)}

   ll_c <- mixll_gamma(a, Gammaest, Sigmaest, Sx)
   awe[i] <- N*ll_c+(log(N)+1.5)*df*2
  }
  return(list(awe=awe, u=which.min(awe)))
}