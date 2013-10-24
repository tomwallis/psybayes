\name{bern_bin}
\alias{bern_bin}
\title{Bin bernoulli trials on a continuous variable x, calculate stats.}
\usage{
  bern_bin(data, breaks, spacing = "equal", x_name = "x",
    y_name = "y", additional_factors = "none", ...)
}
\arguments{
  \item{data}{a data frame containing (at least) a numeric
  vector \code{x} and a binary outcome vector \code{y}.}

  \item{breaks}{a single number giving the number of bins
  desired, or a vector of breakpoints between cells.}

  \item{spacing}{a character vector of either
  \code{"equal"} or \code{"quantile"}, producing equally
  spaced bins along x or bins with approximately equal
  numbers of data points, respectively. Doesn't do anything
  if you pass \code{breaks} as a vector of breakpoints.}

  \item{x_name}{A character vector specifying the name of
  the continuous variable to bin over.}

  \item{y_name}{A character vector specifying the name of
  the binary outcome variable.}

  \item{additional_factors}{A character vector specifying
  any factors in the dataset to split x over.}

  \item{...}{additional arguments passed to
  \link{beta_cis}.}
}
\value{
  a data frame with x bins as rows, and the number of
  successes, fails, and trials in each bin.
}
\description{
  Takes an x and a y vector, where y is binary (bernoulli
  trials).
}
\examples{
library(ggplot2)
dat <- data.frame(x = rnorm(100), y = rbinom(100,1,prob=0.5))
qplot(dat$x,dat$y)
equal_spacing <- bern_bin(dat,breaks=5)
equal_numbers <- bern_bin(dat,breaks=5,spacing='quantile')
custom_spacing <- bern_bin(dat,breaks=seq(-3,3,l=6))

# plot:
equal_spacing$label <- 'Equal spacing'
equal_numbers$label <- 'Equal numbers'
custom_spacing$label <- 'Custom'

plot_dat <- rbind(equal_spacing,equal_numbers,custom_spacing)

fig <- ggplot(plot_dat,aes(x=xmid,xmax=xmax,xmin=xmin,y=y,ymax=ymax,ymin=ymin)) + geom_point() + geom_errorbar() + geom_errorbarh()
fig <- fig + facet_wrap(~ label, ncol=1)
fig

## Show use of additional factors, and changing beta distribution calculations:
dat$a_factor <- rbinom(nrow(dat),1,prob=0.5)
dat$a_factor <- factor(dat$a_factor)

binned <- bern_bin(dat,breaks=5, additional_factors = "a_factor", probs = c(0.025, 0.975))
fig <- ggplot(binned,aes(x=xmid,xmax=xmax,xmin=xmin,y=y,ymax=ymax,ymin=ymin)) + geom_point() + geom_errorbar() + geom_errorbarh()
fig <- fig + facet_wrap(~ a_factor, ncol=1)
fig
}
\author{
  Thomas Wallis
}
\seealso{
  \link{beta_cis}
}

