\name{calc_deviance}
\alias{calc_deviance}
\title{Compute deviance for a binomial trial dataset.}
\usage{
calc_deviance(y, p, size = 1)
}
\arguments{
  \item{y}{vector of successes}

  \item{p}{vector of probabilities (could be predicted
  probabilities from a model)}

  \item{size}{number of trials (1 = bernoulli trial)}
}
\description{
Calculate the model deviance for a vector of successes y
given an equal-length vector of p values and the number of
trials in each row.
}
\details{
The deviance is equal to 2 times the log likelihood of the
model minus the log likelihood of the saturated model (with
a parameter for every data point): \deqn{D = 2 (L_max -
L_model)}

The log likelihood of the saturated model (L_max) is
computed by calculating the log likelihood of a model with
the probabilities set to the value of each data point.

Deviance is a scale where lower numbers are better.
Explaining all the variance (the saturated model) returns a
deviance of 0, whereas the upper bound for deviance is
infinity (observed data are impossible given model
predictions).

Warning: this simple method for calculating the saturated
model deviance may not be appropriate for more complex
models with nesting. See e.g.
\url{http://warnercnr.colostate.edu/~gwhite/mark/markhelp/saturatedmodel.htm}
}
\examples{
y <- rbinom(100,size=1,prob=0.5)
(deviance <- calc_deviance(y=y,p=0.5,size=1))
(saturated_dev <- calc_deviance(y=y,p=y,size=1))
}
\author{
Thomas Wallis
}
\references{
For a discussion of deviance in the context of psychometric
function fitting, see Wichmann, F. A., & Hill, N. J.
(2001). The psychometric function: I. Fitting, sampling,
and goodness of fit. Perception & Psychophysics, 63(8),
1293–1313.
}

