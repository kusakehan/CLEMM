\name{clustering_err}
\alias{clustering_err}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clustering error rate
%%  ~~function to do ... ~~
}
\description{
Calculate clustering error rate by using estimations
}
\usage{
clustering_err(K, da, em_res, pi0, mu, Sigma, idx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{da}{Data for clustering}
  \item{K}{Selected number of clusters}
  \item{em_res}{Results list returned by an EM algorithm}
  \item{pi0}{True mixing weights, if known}
  \item{mu}{True clusters centers, if known}
  \item{Sigma}{True clusters covariance, if known}
  \item{idx}{True clusters labels, if known}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{cluster_err}{clustering error rate}
\item{mean_err}{estimation error of clusters centers}
\item{wt_err}{estimation error of clusters mixing weights}
\item{cov_err}{estimation error of clusters covariance}
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
res_clemm = clemm_em(dat, K, u, iter=800, opts=opts, init=init)

clemm_err = clustering_err(K, dat, em_res=res_clemm, pi0=pi0, mu=da$mu, Sigma=da$Sigma, idx=idx)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
