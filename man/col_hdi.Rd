\name{col_hdi}
\alias{col_hdi}
\title{Calculate the HDI along the columns of a matrix.}
\usage{
  col_hdi(X, cred_mass = 0.95)
}
\arguments{
  \item{X}{a matrix of numbers}

  \item{probs}{vector of quantiles to compute.}
}
\value{
  a matrix where columns correspond to columns of \code{X},
  and rows to the lower and upper HDI limits.
}
\description{
  Receives a matrix X as input, and calculates the HDI
  along the columns using sapply.
}
\examples{
X <- matrix(rnorm(1000),ncol=10)
col_hdi(X)
}
\author{
  Thomas Wallis
}
\seealso{
  Other distribution_summaries: \code{\link{col_quantile}},
  \code{\link{hdi}}
}

