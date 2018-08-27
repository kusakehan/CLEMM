\name{gmm_em}
\alias{gmm_em}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm for GMM
}
\description{
EM algorithm for GMM estimation
}
\usage{
gmm_em(da, K, iter, stopping=1e-7, init, typ="G", add_eye="False")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data for clustering}
  \item{K}{Selected number of clusters}
  \item{iter}{Maximum iteration for EM algorithm}
  \item{stopping}{Stopping criterion}
  \item{init}{Initial values for EM}
  \item{typ}{CLEMM Type, "G": General CLEMM; "S": CLEMM-Shared}
  \item{add_eye}{Add 0.01*I_r to prevent S_k being too singular}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{ll_seq}{Observed Log-likelihood sequence}
\item{a}{Estimation of mixing weights}
\item{muEst}{Estimation of clusters centers}
\item{SigmaEst}{Estimation of clusters covariance}
\item{iteration}{Number of iterations}
%% ...
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
K <- 3;r <- 15; u<- 1;N <- 1000;md=1
da <- data_generation(r, u, N, K, md)
pi0 <- c(0.3, 0.2, 0.5)
idx <- sample(1:K, N, replace = TRUE, prob = pi0)
dat <- matrix(NA, N, r)
for (j in 1:K) {
    if(md <= 3){
      x_tmp <- rmvnorm(sum(idx==j), mean = da$mu[, j], sigma = da$Sigma[, , j])
    }else{ x_tmp <- rmvnorm(sum(idx==j), mean = da$mu[, j], sigma = da$Sigma)}
    dat[idx==j, ] <- x_tmp
}
init=list()
init$centers = da$mu
init$wt = pi0
init$cov = da$Sigma

res_clemm = gmm_em(dat, K, iter=800, init=init)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
