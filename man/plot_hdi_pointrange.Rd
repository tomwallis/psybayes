\name{plot_hdi_pointrange}
\alias{plot_hdi_pointrange}
\title{Mean and HDI of a vector as ggplot2 pointrange.}
\usage{
plot_hdi_pointrange(...)
}
\arguments{
  \item{...}{additional values passed to stat_summary.}
}
\description{
This is a wrapper for ggplot2's stat_summary function that
acts on a y vector from your ggplot2 object, producing a
pointrange geom that you add to a plot. The point shows the
mean of the distribution, the range shows the 95% HDI.
}
\examples{
dat <- data.frame(x = c("Distribution 1", "Distribution 2"),
y = rnorm(100))
dat$y[dat$x == "Distribution 1"] <- dat$y[dat$x == "Distribution 1"] + 1
library(ggplot2)
fig <- ggplot(dat, aes(x = x, y = y)) + plot_hdi_pointrange()
fig
}
\author{
Thomas Wallis
}

