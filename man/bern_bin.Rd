\name{bern_bin}
\alias{bern_bin}
\title{Bin bernoulli trials on a continuous variable x, calculate stats.}
\usage{
  bern_bin(x, y, breaks, equal_spacing = TRUE)
}
\arguments{
  \item{x}{a numeric vector}

  \item{y}{binary outcomes the same length as x}

  \item{breaks}{a single number giving the number of bins
  desired, or a vector of breakpoints between cells.}

  \item{equal_spacing}{if true, bins are equally spaced
  along x. If false, bins are of approximately equal sizes,
  centred on average. Doesn't do anything if you pass
  \code{breaks} as a vector of breakpoints.}
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
x <- rnorm(100)
y <- rbinom(100,1,prob=0.5)
qplot(x,y)
equal_spacing <- bern_bin(x,y,breaks=5)
equal_numbers <- bern_bin(x,y,breaks=5,equal_spacing=FALSE)
custom_spacing <- bern_bin(x,y,breaks=seq(-3,3,l=6))

# beta distribution confidence intervals:
equal_spacing <- beta_cis(equal_spacing)
equal_numbers <- beta_cis(equal_numbers)
custom_spacing <- beta_cis(custom_spacing)

# plot:
equal_spacing$label <- 'Equal spacing'
equal_numbers$label <- 'Equal numbers'
custom_spacing$label <- 'Custom'

plot_dat <- rbind(equal_spacing,equal_numbers,custom_spacing)

fig <- ggplot(plot_dat,aes(x=xmid,xmax=xmax,xmin=xmin,y=y,ymax=ymax,ymin=ymin)) + geom_point() + geom_errorbar() + geom_errorbarh()
fig <- fig + facet_wrap(~ label, ncol=1)
fig
}
\author{
  Thomas Wallis
}

