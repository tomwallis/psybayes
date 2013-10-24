\name{beta_cis}
\alias{beta_cis}
\title{From binomial data frame, compute confidence limits of beta distribution.}
\usage{
  beta_cis(data, probs = c(0.025, 0.5, 0.975),
    rule_of_succession = TRUE)
}
\arguments{
  \item{data}{a data frame.}

  \item{probs}{probabilities to return (increasing). If
  this is a length 3 vector, the values will be assigned to
  ymin, y and ymax. If it's a vector of length 2, then this
  specifies ymin and ymax, and y will be the proportion
  correct.}

  \item{rule_of_succession}{If true, add one success and
  one failure to every observation (Laplace's rule of
  succession correction).}
}
\value{
  the data frame with added y, ymin and ymax columns. Data
  (n_successes, etc) will be original, not
  rule-of-succession corrected.
}
\description{
  Input is a data frame containing columns n_success,
  n_fails.
}
\examples{
example(bern_bin)
}
\author{
  Thomas Wallis
}
\seealso{
  \link{bern_bin}
}

