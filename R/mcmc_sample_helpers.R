

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