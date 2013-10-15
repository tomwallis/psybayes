\name{hdi}
\alias{hdi}
\title{Calculate the HDI (highest density interval) of a vector.}
\usage{
  hdi(x, cred_mass = 0.95)
}
\arguments{
  \item{x}{a vector of numbers.}

  \item{cred_mass}{the credible mass to calculate (scalar
  from 0 to 1).}
}
\value{
  a vector of length 2 containing the limits of the HDI.
}
\description{
  Receives a vector x as input, and calculates the hdi for
  a given mass.
}
\examples{
x <- rnorm(1000)
hdi(x)
}
\author{
  Thomas Wallis modified original code by John Kruschke,
  from Doing Bayesian Data Analysis.
}
\seealso{
  \link{col_quantile} \link{col_hdi}
}

