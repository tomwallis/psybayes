## functions for calculating binomial deviance and other metrics given y (true binomial result) and yhat (predicted probability) inputs.

# TSAW. tsawallis@gmail.com

# Log likelihood ---------------------------------------------------
#' Log likelihood for a bernoulli trial dataset.
#' 
#' This function calculates the log likelihood for a vector of binomial observations y
#' given an equal-length vector of p values and the number of trials in each row.
#' This can be useful, for example, when computing log likelihood on a crossvalidated dataset,
#' where y is the observed responses and p are the model predictions.
#' 
#' @export  
#' 
#' @param y vector of successes
#' @param p vector of probabilities
#' @param size number of trials (1 = bernoulli trial)
#' @author Thomas Wallis
#' @examples
#' y <- rbinom(n=10,size=1,prob=.6)
#' yhat <- c(0.4,0.5,0.5,0.5,.6,.6,.6,.6,.7,.7)
#' (ll <- calc_ll(y,yhat,size=1))

calc_ll <- function(y,p,size=1){
  # Manually:
#   log_like <- lchoose(size,y) + y * log(p) + (size - y) * log(1-p)
#   # sum over logs to get the final likelihood:
#   total_log_like <- sum(log_like)  
#   return(total_log_like)
  
  # using R density:
  return(sum(dbinom(y,size=size,prob=p,log=TRUE)))
}


# Deviance ---------------------------------------------------
#' Compute deviance for a binomial trial dataset.
#' 
#' Calculate the model deviance for a vector of successes y
#' given an equal-length vector of p values and the number of trials in each row.
#' 
#' @export  
#' 
#' @param y vector of successes
#' @param p vector of probabilities (could be predicted probabilities from a model)
#' @param size number of trials (1 = bernoulli trial)
#' @author Thomas Wallis
#' 
#' @details The deviance is equal to -2 times the log likelihood of the model minus the 
#' log likelihood of the saturated model (with a parameter for every data point):
#' \deqn{D = -2 (L_model - L_saturated}
#' 
#' The log likelihood of the saturated model is computed by calculating the log likelihood
#' of a model with the probabilities set to the value of each data point.
#' 
#' Deviance is a scale where lower numbers are better. Explaining all the variance (the
#' saturated model) returns a deviance of 0.
#' 
#' Warning: this simple method for calculating the saturated model deviance may not be
#' appropriate for more complex models with nesting. See e.g. \url{http://warnercnr.colostate.edu/~gwhite/mark/markhelp/saturatedmodel.htm}
#' 
#' @examples
#' y <- rbinom(100,size=1,prob=0.5)
#' (deviance <- calc_deviance(y=y,p=0.5,size=1))
#' (saturated_dev <- calc_deviance(y=y,p=y,size=1))

calc_deviance <- function(y,p,size=1){
  # saturated model has one param per data point:
  saturated_ll <- calc_ll(y=y,p=(y/size),size=size)
  
  model_ll <- calc_ll(y=y,p=p,size=size)

  # deviance = -2 times the log likelihood differences.  
  dev <- -2 * (model_ll - saturated_ll)
  return(dev)
}


# WAIC ---------------------------------------------------
# as defined in Gelman, Hwang & Vehtari, 2013. Understanding predictive information criteria for Bayesian models.  (Equations 12 -- 13).




