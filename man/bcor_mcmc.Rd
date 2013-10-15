\name{bcor_mcmc}
\alias{bcor_mcmc}
\title{Estimate Bayesian (posterior) correlation coefficients.}
\usage{
  bcor_mcmc(data_frame, prior = "lkj_corr(1.0)", ...)
}
\arguments{
  \item{data_frame}{A data frame where rows are
  observations and columns are the variables whose
  correlation structure you wish to estimate. All variables
  must be continuous (no factors).}

  \item{prior}{A text string containing the (Stan model)
  prior specification for the correlation matrices.}

  \item{...}{Additional variables passed to stan_sample
  function.}
}
\value{
  Returns an rstan fit object.
}
\description{
  Estimate Pearson's r and Spearman's rho between the
  columns in a data frame. This is a wrapper function for
  Stan.
}
\details{
  This function takes a numerical data matrix X as input.
  From this, the Pearson r and Spearman's rho (rank order
  correlation coefficient) will be estimated using MCMC
  sampling from the \link{stan} package. Currently the
  function standardises the data internally (using the
  \link{scale} function), but could be easily modified to
  also estimate the mean and variance of each variable.

  The \code{prior} is applied to each correlation matrix
  (Pearson and Spearman). The default is a multivariate
  version of a Beta distribution with parameters (1,1),
  i.e. it is multivariate uniform, scaled over the range -1
  to 1. See Stan manual for details.

  For a big example of all this in action, run
  \code{example(bcor_mcmc)}.
}
\examples{
# generate random variables with a specified correlation structure (from: http://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/)
# n_var variables. Assume that they are standardised (mean 0 sd 1).
n_vars <- 4
# number of observations to simulate
n_obs = 200
# correlation matrix:
M = matrix(c(1, 0.7, 0.7, 0.5,
             0.7, 1, 0.95, 0.3,
             0.7, 0.95, 1, 0.3,
             0.5, 0.3, 0.3, 1), nrow=4, ncol=4)

# Cholesky decomposition
L <- chol(M)
nvars = dim(L)[1]
# Random variables that follow an M correlation matrix
r <- t(L) \%*\% matrix(rnorm(n_vars*n_obs), nrow=n_vars, ncol=n_obs)
r <- t(r)
dat <- data.frame(r)

# Now that we've generated data, fit and analyse:
fit <- bcor_mcmc(dat)
print(fit)
bcor_plot(dat,fit,true_correlation=M)
tables <- bcor_table(dat,fit)
tables$pearson_table
cor(dat)
}
\author{
  Thomas Wallis.
}
\seealso{
  \link{bcor_plot}, \link{bcor_table}, \link{stan_sample},
  \link{stan}
}

