\name{bcor_table}
\alias{bcor_table}
\title{Table of correlation coefficients with HDIs.}
\usage{
  bcor_table(data_frame, fit, variable_names = FALSE,
    central_tendency_function = mean,
    spread_function = hdi)
}
\arguments{
  \item{data_frame}{A data frame where rows are
  observations and columns are the variables whose
  correlation structure you wish to estimate. All variables
  must be continuous (no factors).}

  \item{fit}{The stanfit object produced by
  \link{bcor_mcmc}.}

  \item{variable_names}{If true, include variable names
  from the data_frame in table.}

  \item{central_tendency_function}{A function to calculate
  the data's central tendency. Should return a scalar with
  input as a vector.}

  \item{spread_function}{A function to calculate the data's
  spread (e.g. hdi or quantiles). Should return two values
  (lower, upper) given input as a vector.}
}
\value{
  A list of tables suitable for passing to the xtable
  function. One table is for the Pearson Product Moment
  correlation coefficent, one is for Spearman's rho (rank
  order correlation) and the third are the same computed
  from samples of the prior distribution.
}
\description{
  Generate a numerical matrix of the correlation
  coefficients estimated by bcor_mcmc.
}
\examples{
See the example for bcor_mcmc.
}
\author{
  Thomas Wallis.
}
\seealso{
  \link{bcor_mcmc}, \link{bcor_plot}
}

