\name{col_quantile}
\alias{col_quantile}
\title{Calculate quantiles along the columns of a matrix.}
\usage{
  col_quantile(X, probs = c(0.025, 0.975))
}
\arguments{
  \item{X}{a matrix of numbers}

  \item{probs}{vector of quantiles to compute.}
}
\value{
  a matrix where columns correspond to columns of \code{X},
  and rows to each requested \code{prob}.
}
\description{
  Receives a matrix X as input, and calculates quantiles
  along the columns using sapply.
}
\examples{
X <- matrix(rnorm(1000),ncol=10)
col_quantile(X)
}
\author{
  Thomas Wallis
}
\seealso{
  \link{hdi}
}

