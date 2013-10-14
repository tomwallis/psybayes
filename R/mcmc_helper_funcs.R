# A bunch of generic helper functions for analysing mcmc output.
# TSAW wrote it, unless otherwise noted. tsawallis@gmail.com

# Vectorised quantile function-------------------
vectorQuants <- function(x,probs=c(0.025,0.975)){
  n_cols <- NCOL(x)
  q <- matrix(rep(NA,times=n_cols*length(probs)),nrow=length(probs))
  for (i in 1:n_cols){
    q[,i] <- quantile(x[,i],probs=probs)
  }
  return(q)
}

## HDI of mcmc---------------------------

hdiOfMCMC = function( sampleVec , credMass=0.95 ) {
  # Function by John Kruschke, from "Doing Bayesian Data Analysis" text.
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

# vectorised hdi---------------------
# 
vectorHDI <- function(x,credMass=0.95){
  n_cols <- NCOL(x)
  q <- matrix(rep(NA,times=n_cols*2),nrow=2)
  for (i in 1:n_cols){
    q[,i] <- hdiOfMCMC(x[,i],credMass=credMass)
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