# A bunch of generic helper functions for analysing mcmc output.
# TSAW wrote it, unless otherwise noted. tsawallis@gmail.com


# Column-wise quantile function ---------------------------------------------------
#' Calculate quantiles along the columns of a matrix.
#' 
#' Receives a matrix X as input, and calculates quantiles along the columns using sapply.
#' 
#' @export  
#' 
#' @param X a matrix of numbers
#' @param probs vector of quantiles to compute.
#' @return a matrix where columns correspond to columns of \code{X}, and rows to each requested \code{prob}.
#' @author Thomas Wallis
#' @seealso \link{hdi}
#' @examples
#' X <- matrix(rnorm(1000),ncol=10)
#' col_quantile(X)

col_quantile <- function(X,probs=c(0.025,0.975)){
  n_cols <- NCOL(X)
  fun <- function(i,X,probs) quantile(X[,i],probs=probs)
  q <- sapply(1:n_cols,fun,X,probs=probs)
  return(q)
}

# HDI of a vector---------------------------
#' Calculate the HDI (highest density interval) of a vector.
#' 
#' Receives a vector x as input, and calculates the hdi for a given mass.
#' 
#' @export  
#' 
#' @param x a vector of numbers.
#' @param cred_mass the credible mass to calculate (scalar from 0 to 1).
#' @return a vector of length 2 containing the limits of the HDI.
#' @author Thomas Wallis modified original code by John Kruschke, from Doing Bayesian Data Analysis.
#' @seealso \link{col_quantile}
#' @examples
#' x <- rnorm(1000)
#' hdi(x)
#' 
hdi = function( x , cred_mass=0.95 ) {
  sorted_x = sort( x )
  ci_ind = floor( cred_mass * length( sorted_x ) )
  n_ci = length( sorted_x ) - ci_ind
  ci_width = rep( 0 , n_ci )
  
  # iterate through each lower bound, finding the smallest distance.
  fun <- function(i,x,ci_ind) x[i + ci_ind] - x[i]
  ci_width <- sapply(1:n_ci, fun, sorted_x, ci_ind)

  hdi_min = sorted_x[ which.min( ci_width ) ]
  hdi_max = sorted_x[ which.min( ci_width ) + ci_ind ]
  hdi_lim = c( hdi_min , hdi_max )
  return( hdi_lim )
}

# vectorised hdi---------------------
# 
vectorHDI <- function(x,cred_mass=0.95){
  n_cols <- NCOL(x)
  q <- matrix(rep(NA,times=n_cols*2),nrow=2)
  for (i in 1:n_cols){
    q[,i] <- hdiOfMCMC(x[,i],cred_mass=cred_mass)
  }
  return(q)
}



# Wrapper function for fitting stan models---------------------------

stan_sample <- function(model_string,stan_data,
                        n_saved_samples=1000,
                        iter=1000,
                        warmup=iter/2,
                        num_chains=4,
                        thin=floor((iter - warmup)/ (n_saved_samples / num_chains)),
                        seed=sample.int(.Machine$integer.max, 1),
                        parallel=TRUE) {
  
  # Returns a fit object from the stan sampling with the above parameters----------------------.
  
  require(rstan)
  
  # first fit bug test run:
  initial_fit <- stan(model_code = model_string, data = stan_data,iter = 10, chains = 4, seed=seed)
  
  
  if(parallel==TRUE){
    # set up for parallel sampling:
    require(doMC)
    options(cores=4) #set this appropriate to your system
    registerDoMC()
    
    #prep for parallel case
    fun <- function(i){
      fit <- stan(fit = initial_fit, data = stan_data, iter = iter, warmup = warmup, chains = 1, verbose = F, thin=thin, chain_id=i, seed=seed)
      return(fit)
    }
    
    #parallel case
    start <- proc.time()[3]
    parallel_fit <- foreach(i = 1:num_chains) %dopar% fun(i)
    (proc.time()[3] - start)/60
    
    # squish fit objects together:
    fit <- sflist2stanfit(parallel_fit)
    
  } else {
    fit <- stan(fit = initial_fit, data = stan_data, iter = iter, warmup = warmup, chains = num_chains, verbose = F, thin=thin, seed=seed)
    
  }
  
  return(fit) 
  
}