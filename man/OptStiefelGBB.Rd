\name{OptStiefelGBB}
\alias{OptStiefelGBB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A feasible method for optimization with orthogonality constraints full
}
\description{
Curvilinear search algorithm for optimization on Stiefel manifold
}
\usage{
OptManiMulitBallGBB(X, opts, fun, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{n by k matrix such that X'*X = I}
  \item{fun}{objective function and its gradient:
  #      [F, G] = fun(X,  data1, data2)
  #      F, G are the objective function value and gradient, repectively
  #      data1, data2 are addtional data, and can be more}
\item{opts}{option structure with fields:
  #             record = 0, no print out
  #             mxitr       max number of iterations
  #             xtol        stop control for ||X_k - X_{k-1}||
  #             gtol        stop control for the projected gradient
  #             ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
  #                         usually, max{xtol, gtol} > ftol}
}
\details{
Used in clemm_em function
}
\value{
 \item{X}{solution}
 \item{Out}{output information}
}
\references{
Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1-2), 397-434.
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
Gammaest <- OptStiefelGBB(Gammaest_init, opts, FGfun, p, Sx, S_tmp)$X
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
