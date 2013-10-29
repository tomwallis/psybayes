\name{bern2binom}
\alias{bern2binom}
\title{Binomial design matrix from Bernoulli trials.}
\usage{
  bern2binom(data, formula, group_factors = "none",
    response = "y")
}
\arguments{
  \item{data}{A data frame with rows as Bernoulli trials.}

  \item{formula}{The formula for the model design at the
  lowest level of the model.}

  \item{group_factor}{A character vector specifying the
  name of one or more grouping variables found in
  \code{data}. Usually "subject".}

  \item{response}{A character vector specifying the name of
  the response variable (correct or incorrect).}
}
\value{
  Returns a new data frame.
}
\description{
  Transform a data frame where each row is a Bernoulli
  trial into a binomial data frame (ntrials, nsucceses)
  based on the design matrix for a model.
}
\details{
  In Stan it is often much more efficient to sample a model
  binomially rather than as individual Bernoulli trials.
  This function can be used to transform a data frame of
  Bernoulli trials (like you might save in an experiment)
  into a binomial data frame, based on a formula for the
  model.
}
\examples{
data(detection)
head(detection)
# detection has rows as Bernoulli trials.
# we would like to model performance as a function of log contrast, for each subject
# within a population.
binom_list <- bern2binom(detection, correct ~ log(contrast), group_factors = "subject", response = "correct")

binom_list$binom_dat

# design matrix can be passed to Stan along with response and grouping variables:
binom_list$design_matrix
}
\author{
  Thomas Wallis.
}
\seealso{
  Other data_munging: \code{\link{bern_bin}}
}

