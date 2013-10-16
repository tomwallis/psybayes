\name{bcor_plot}
\alias{bcor_plot}
\title{Plot posterior correlation coefficients.}
\usage{
  bcor_plot(data_frame, fit, out_file = "none",
    true_correlations = FALSE, plot_prior = TRUE)
}
\arguments{
  \item{data_frame}{A data frame where rows are
  observations and columns are the variables whose
  correlation structure you wish to estimate. All variables
  must be continuous (no factors).}

  \item{fit}{The stanfit object produced by
  \link{bcor_mcmc}.}

  \item{out_file}{The filename of the output file (with
  .pdf extension), or 'none' (print to R directly).}

  \item{true_correlations}{If a matrix of numbers is passed
  here, this is taken to be the True correlation
  coefficients used to generate the data. If left at FALSE,
  the generating correlation coefficients are not plotted.}

  \item{plot_prior}{If TRUE, plot an estimate of the prior
  density from the model.}
}
\value{
  If \code{out_file} contains a filename with .pdf
  extension, a plot of each pairwise correlation is sent to
  out_file. Otherwise these are printed to R. The dark
  histogram shows samples from the posterior, the faint
  histogram shows samples from the prior.
}
\description{
  Plot the correlation coefficients estimated by bcor_mcmc.
}
\examples{
See the example for bcor_mcmc.
}
\author{
  Thomas Wallis.
}
\seealso{
  Other bcor: \code{\link{bcor_mcmc}},
  \code{\link{bcor_table}}
}

