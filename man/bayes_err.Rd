\name{bayes_err}
\alias{bayes_err}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayes error
}
\description{
Calculate Bayes error of using true parameters
}
\usage{
bayes_err(da, pi0, mu, Sigma, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{da}{Data for clustering}
  \item{pi0}{True clusters weights}
  \item{mu}{True clusters centers}
  \item{sigma}{True clusters covariance}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Bayes' classification error
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

## First generate data
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

bayeserr(dat, pi0, da$mu, da$Sigma, idx)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
