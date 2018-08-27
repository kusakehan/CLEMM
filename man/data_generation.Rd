\name{data_generation}
\alias{data_generation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data generation
%%  ~~function to do ... ~~
}
\description{
Generate model (M1)-(M5) from the simulation Section of the paper
}
\usage{
data_generation(r, u, N, K, md, seed=999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{Dimension of data}
  \item{u}{Envelope dimension}
  \item{K}{Selected number of clusters}
  \item{md}{Type of model in the paper, choose from 1 to 5}
  \item{seed}{Random seed for generation process, default 999}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{mu}{Clusters centers}
\item{Sigma}{Clusters covariance}
\item{Gamma}{Envelope subspace basis}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
