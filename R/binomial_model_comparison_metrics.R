## functions for calculating binomial deviance and other metrics given y (true binomial result) and yhat (predicted probability) inputs.

# TSAW. tsawallis@gmail.com

# Log likelihood ---------------------------------------------------
#' Log likelihood for a bernoulli trial dataset.
#' 
#' This function calculates the log likelihood for a vector of binomial observations
#' given an equal-length vector of p values and the number of trials in each row.
#' 
#' @export  
#' 
#' @param y vector of successes
#' @param p vector of probabilities
#' @param n vector of trial numbers (1 = bernoulli trial)
#' @author Thomas Wallis
#' @examples
#' y <- rbinom(n=10,size=1,prob=.6)
#' yhat <- c(0.4,0.5,0.5,0.5,.6,.6,.6,.6,.7,.7)
#' (ll <- calc_ll(y,yhat,n=1))

calc_ll <- function(y,p,n=1){
  # Manually:
#   log_like <- lchoose(n,y) + y * log(p) + (n - y) * log(1-p)
#   # sum over logs to get the final likelihood:
#   total_log_like <- sum(log_like)  
#   return(total_log_like)
  
  # using R density:
  return(sum(dbinom(y,size=n,prob=p,log=TRUE)))
}


# Deviance ---------------------------------------------------
#' Compute deviance for a bernoulli trial dataset.
#' 
#' This function calculates the log likelihood for a vector of binomial observations
#' given an equal-length vector of p values and the number of trials in each row.
#' 
#' @export  
#' 
#' @param y vector of successes
#' @param p vector of probabilities
#' @param n vector of trial numbers (1 = bernoulli trial)
#' @author Thomas Wallis
#' @examples
#' y <- rbinom(n=10,size=1,prob=.6)
#' yhat <- c(0.4,0.5,0.5,0.5,.6,.6,.6,.6,.7,.7)
#' (ll <- calc_ll(y,yhat,n=1))
calc_deviance <- function(y,p,n=1){
  # log like * -2:
  log_like <- calc_ll(y=y,p=p,n=n)
  return(log_like*-2)
}


# WAIC ---------------------------------------------------
# as defined in Gelman, Hwang & Vehtari, 2013. Understanding predictive information criteria for Bayesian models.  (Equations 12 -- 13).




