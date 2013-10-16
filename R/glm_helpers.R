# Linear predictor from a matrix of betas ----------------------------------
#' Calculate the linear predictor.
#' 
#' Given a design matrix X and a matrix of betas, calculate the linear 
#' predictor for each row-wise entry of beta. If beta is a vector, returns the 
#' single vector linear predictor.
#' 
#' @export  
#' 
#' @param X       Design matrix. Should have \code{ncol(X) == ncol(beta)}.
#' @param beta   A matrix of credible beta weights (columns = beta_1 ... beta_n,
#' rows = samples from MCMC), or a vector of beta weights.
#' 
#' @details
#' In mcmc applications, we don't end up with a point estimate of beta weights
#' from a GLM but a vector of credible samples for each beta weight. To compute
#' credible curves on the response scale (i.e. after running through the link function),
#' we must first compute the values of the linear predictor for each MCMC sample.
#' 
#' This function does so using sapply.
#' 
#' Alternatively, if \code{beta} is a vector, the function simply returns the 
#' inner product of X and beta.
#' 
#' @return Returns a matrix with rows = samples and columns = y(x). I've arranged it in this way
#' to make it easier to examine means of y(x), e.g. with \link{colMeans}, and calculate
#' hdis (with \link{hdi}). 
#' 
#' Alternatively, a vector of the inner product of X and beta.
#' 
#' @family glm_helpers
#' 
#' @author Thomas Wallis.
#' @examples
#' X <- matrix(rnorm(1000),ncol=10)
#' betas <- matrix(rep(1:200),ncol=10) # like we've taken 20 mcmc samples.
#' eta <- linpred(X,betas)

linpred <- function(X, beta){
  if(is.matrix(beta)){
    if(NCOL(X) != NCOL(beta)) stop("columns on X and beta do not match")  
    fun <- function(i,X,beta) X %*% beta[i,]
    eta <- sapply(1:nrow(beta),fun,X,beta)
    return(t(eta))    
  } else {
    eta <- X %*% beta
    eta <- as.vector(eta)
    return(eta)
  }
}



# Define link functions ----------------------------------
#' Define link functions.
#' 
#' This function can be used to create an S3 object of class "psybayes_link" that
#' includes link and inverse link functions for transforming from linear predictor
#' space into probability space, as well as code for inserting into Stan models.
#' It is a modified version of \link{make.link} from the stats package.
#' 
#' @export
#' @import wutils ggplot2  
#' 
#' @param link   Character; one of "logit", "probit", "cauchy", "cloglog".
#' 
#' @return Returns an object of class "psybayes_link". This is really just a list, with the 
#' following components:
#' \itemize{
#' \item \code{linkfun}
#' \item \code{linkinv}
#' \item \code{stan_string}
#' \item \code{name}
#' }
#' 
#' @details
#' The link functions returned in the object accept three arguments:
#' \itemize{
#' \item \code{mu} or \code{eta}: for link and inverse link respectively; the p value or the linear predictor.
#' \item \code{gamma}: the lower bound of the psychometric function (defaults to zero).
#' \item \code{lambda}: the lapse rate of the psychometric function (upper bound is \code{1 - lambda}; defaults to zero).
#' }
#' 
#' @author Thomas Wallis.
#' @family glm_helpers
#' @examples
#' logistic <- psy_link('logit')
#' # find logit for 0.5:
#' logistic$linkfun(0.5)
#' 
#' library(ggplot2)
#' library(wutils)
#' 
#' # examine link functions:
#' links <- c('logit','probit','cauchy','cloglog')
#' 
#' # there is probably a nicer way to do this with plyr but I can't figure it right now:
#' d1 <- data.frame()
#' d2 <- data.frame()
#' for(i in 1:length(links)){
#'   link_fun <- psy_link(links[i])
#'   
#'   # link function:
#'   this_dat <- edply(list(mu=seq(0.01,0.99,l=200),gamma=c(0,0.25,0.5),lambda=c(0,0.05)),link_fun$linkfun)
#'   this_dat$link <- links[i]
#'   d1 <- rbind(d1,this_dat)
#'   
#'   # inverse link:
#'   this_dat <- edply(list(eta=seq(-10,10,l=200),gamma=c(0,0.25,0.5),lambda=c(0,0.05)),link_fun$linkinv)
#'   this_dat$link <- links[i]
#'   d2 <- rbind(d2,this_dat)
#' }
#' 
#' ggplot(d1,aes(x=mu,y=V1,colour=link)) + geom_line() + facet_grid(gamma ~ lambda) + xlab('p(c)') + ylab('linear predictor')
#' ggplot(d2,aes(x=eta,y=V1,colour=link)) + geom_line() + facet_grid(gamma ~ lambda) + xlab('linear predictor') + ylab('p(c)')


psy_link <- function(link)
{
  switch(link, 
         logit = {
           linkfun <- function(mu, gamma=0, lambda=0) {
             p <- (gamma - mu) / (gamma + lambda - 1)
             return(qlogis(p))
           }
           linkinv <- function(eta, gamma=0, lambda=0) {
             gamma + (1 - gamma - lambda) * plogis(eta)
           }
         }, 
         probit = {
           linkfun <- function(mu, gamma=0, lambda=0) {
             p <- (gamma - mu) / (gamma + lambda - 1)
             return(qnorm(p))
           }
           linkinv <- function(eta, gamma=0, lambda=0) {
             thresh <- -qnorm(.Machine$double.eps)
             eta <- pmin(pmax(eta, -thresh), thresh)
             return( gamma + (1 - gamma - lambda) * pnorm(eta))
           }
         },
         cauchy = {
           linkfun <- function(mu, gamma=0, lambda=0) {
             p <- (gamma - mu) / (gamma + lambda - 1)
             qcauchy(p)
           }
             
           linkinv <- function(eta, gamma=0, lambda=0) {
             thresh <- -qcauchy(.Machine$double.eps)
             eta <- pmin(pmax(eta, -thresh), thresh)
             return( gamma + (1 - gamma - lambda) * pcauchy(eta))
           }
         }, 
         cloglog = {
           linkfun <- function(mu, gamma=0, lambda=0) {
             p <- (gamma - mu) / (gamma + lambda - 1)
             log(-log(1 - p))
           } 
           linkinv <- function(eta, gamma=0, lambda=0) {
             eta <- pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
             return( gamma + (1 - gamma - lambda) * eta)
           } 
         }, 
#          weibull = {
#            linkfun <- function(mu, gamma=0, lambda=0) {
#              p <- (gamma - mu) / (gamma + lambda - 1)
#              return(qweibull(p))
#            }
#            linkinv <- function(eta, gamma=0, lambda=0) {
#              gamma + (1 - gamma - lambda) * pweibull(eta)
#            }
#          }, 
         stop(gettextf("%s link not recognised", sQuote(link)), 
                 domain = NA))
  structure(list(linkfun = linkfun, linkinv = linkinv, name = link), class = "psybayes_link")
}



