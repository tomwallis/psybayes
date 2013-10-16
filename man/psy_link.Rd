\name{psy_link}
\alias{psy_link}
\title{Define link functions.}
\usage{
  psy_link(link)
}
\arguments{
  \item{link}{Character; one of "logit", "probit",
  "cauchy", "cloglog".}
}
\value{
  Returns an object of class "psybayes_link". This is
  really just a list, with the following components:
  \itemize{ \item \code{linkfun} \item \code{linkinv} \item
  \code{stan_string} \item \code{name} }
}
\description{
  This function can be used to create an S3 object of class
  "psybayes_link" that includes link and inverse link
  functions for transforming from linear predictor space
  into probability space, as well as code for inserting
  into Stan models. It is a modified version of
  \link{make.link} from the stats package.
}
\details{
  The link functions returned in the object accept three
  arguments: \itemize{ \item \code{mu} or \code{eta}: for
  link and inverse link respectively; the p value or the
  linear predictor. \item \code{gamma}: the lower bound of
  the psychometric function (defaults to zero). \item
  \code{lambda}: the lapse rate of the psychometric
  function (upper bound is \code{1 - lambda}; defaults to
  zero). }
}
\examples{
logistic <- psy_link('logit')
# find logit for 0.5:
logistic$linkfun(0.5)

library(ggplot2)
library(wutils)

# examine link functions:
links <- c('logit','probit','cauchy','cloglog')

# there is probably a nicer way to do this with plyr but I can't figure it right now:
d1 <- data.frame()
d2 <- data.frame()
for(i in 1:length(links)){
  link_fun <- psy_link(links[i])

  # link function:
  this_dat <- edply(list(mu=seq(0.01,0.99,l=200),gamma=c(0,0.25,0.5),lambda=c(0,0.05)),link_fun$linkfun)
  this_dat$link <- links[i]
  d1 <- rbind(d1,this_dat)

  # inverse link:
  this_dat <- edply(list(eta=seq(-10,10,l=200),gamma=c(0,0.25,0.5),lambda=c(0,0.05)),link_fun$linkinv)
  this_dat$link <- links[i]
  d2 <- rbind(d2,this_dat)
}

ggplot(d1,aes(x=mu,y=V1,colour=link)) + geom_line() + facet_grid(gamma ~ lambda) + xlab('p(c)') + ylab('linear predictor')
ggplot(d2,aes(x=eta,y=V1,colour=link)) + geom_line() + facet_grid(gamma ~ lambda) + xlab('linear predictor') + ylab('p(c)')
}
\author{
  Thomas Wallis.
}
\seealso{
  Other glm_helpers: \code{\link{linpred}}
}

