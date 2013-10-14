\name{calc_ll}
\alias{calc_ll}
\title{Log likelihood for a bernoulli trial dataset.}
\usage{
  calc_ll(y, p, n = 1)
}
\arguments{
  \item{y}{vector of successes}

  \item{p}{vector of probabilities}

  \item{n}{vector of trial numbers (1 = bernoulli trial)}
}
\description{
  This function calculates the log likelihood for a vector
  of binomial observations given an equal-length vector of
  p values and the number of trials in each row.
}
\examples{
y <- rbinom(n=10,size=1,prob=.6)
yhat <- c(0.4,0.5,0.5,0.5,.6,.6,.6,.6,.7,.7)
(ll <- calc_ll(y,yhat,n=1))
}
\author{
  Thomas Wallis
}

