\name{env_dim_selection}
\alias{env_dim_selection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Envelope dimension selection
}
\description{
Select envelope dimension from a given range
}
\usage{
env_dim_selection(dim_rng, da, K, iter, stopping=1e-7, opts=NULL, init, typ="G", add_eye="False")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dim_rng}{The range of envelope dimension we want to choose from}
  \item{da}{Data for clustering}
  \item{K}{Selected number of clusters}
  \item{iter}{Maximum iteration for EM algorithm}
  \item{stopping}{Stopping criterion}
  \item{init}{Initial values for EM}
  \item{typ}{CLEMM Type, "G": General CLEMM; "S": CLEMM-Shared}
  \item{add_eye}{Add 0.01*I_r to prevent S_k being too singular}
  \item{opts}{Parameters for mannifold optimization}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{awe}{AWE values for different envelope dimensions}
\item{u}{Selected envelope dimension which corresponds to the minimum AWE values}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(999)
K <- 4;r <- 15; u<- 2;N <- 1000;md=2
da <- data_generation(r, u, N, K, md)
pi0 <- rep(0.25, 4)

set.seed(1)
idx <- sample(1:K, N, replace = TRUE, prob = pi0)
dat <- matrix(NA, N, r)
for (j in 1:K) {
    if(md <= 3){
      x_tmp <- rmvnorm(sum(idx==j), mean = da$mu[, j], sigma = da$Sigma[, , j])
    }else{ x_tmp <- rmvnorm(sum(idx==j), mean = da$mu[, j], sigma = da$Sigma)}
    dat[idx==j, ] <- x_tmp
}

opts <- c()
opts$record <- 0
opts$mxitr <- 2000
opts$xtol <- 1e-10
opts$gtol <- 1e-10
opts$ftol <- 1e-10

init=list()
init$centers = da$mu
init$wt = pi0
init$cov = da$Sigma

dim_sel_G = env_dim_selection(1:12, dat, K, iter=800, stopping=1e-7, init=init)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
