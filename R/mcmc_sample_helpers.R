
# stan_sample ---------------------------
#' Wrapper function for fitting stan models. 
#' 
#' A wrapper function for fitting Stan models using Tom's common settings.
#' Also includes some ability to paralellise code. 
#' Requires the rstan package to run.
#' 
#' @export  
#' @import rstan doMC
#' 
#' @param file                A character string file name or connection that R supports containing the text of a model specification in the Stan language.
#' @param model_string        A character vector of a stan model string (beginning with data block, etc.), if file not provided.
#' @param stan_data           A list object containing named data that the stan model expects.
#' @param n_saved_samples     The number of samples to save, across all chains.
#' @param iter                The number of iterations to run the sampler per chain (must be at least n_saved_samples/n_chains).
#' @param warmup              The number of iterations that should be warmup.
#' @param num_chains          The number of independent MCMC chains to run.
#' @param thin                Stan will save every \code{thin} sample.
#' @param seed                The integer seed to use for seeding the chains.
#' @param parallel            Whether to run the model in parallel (one chain per core; uses doMC package).
#' @param n_cores             The number of cores to use on the machine you're running. 
#' @param ...                 Additional arguments passed to \code{stan} function.
#'
#' @return Returns an rstan fit object.
#' @author Thomas Wallis.
#' @seealso \link{stan}
#' @examples
#' data(cars)
#' model <- '
#' data{
#' int N;
#' vector[N] dist;
#' vector[N] speed;
#' }
#' 
#' parameters{
#' vector[2] beta;
#' real<lower=0> sigma;
#' }
#' 
#' model{
#' vector[N] mu;
#' beta ~ normal(0,100);
#' sigma ~ gamma(2,1e-5);
#' mu <- beta[1] + beta[2] * speed;
#' dist ~ normal(mu, sigma);
#' }
#' '
#' 
#' data <- list(
#'   N = nrow(cars) ,
#'   dist = cars$dist ,
#'   speed = cars$speed
#' )
#' 
#' # fit a linear model predicting distance from speed
#' fit <- stan_sample(model_string=model, stan_data=data)
#' print(fit)
#' 
#' # compare to lm:
#' summary(lm(dist ~ speed, cars))

stan_sample <- function(file,
                        model_string = "",
                        stan_data = list(),
                        n_saved_samples=1000,
                        iter=1000,
                        warmup=floor(iter/2),
                        num_chains=4,
                        thin=floor((iter - warmup)/ (n_saved_samples / num_chains)),
                        seed=sample.int(.Machine$integer.max, 1),
                        parallel=TRUE,
                        n_cores=4,
                        ...) {
  
  # Returns a fit object from the stan sampling with the above parameters----------------------.
  
  require(rstan)
  
  # if supplied as model string, else as file:
  if (model_string == ""){
    # first fit bug test run:
    initial_fit <- stan(file = file, data = stan_data,iter = 10, chains = 4, seed=seed)
  }
  
  if (model_string != ""){
    # first fit bug test run:
    initial_fit <- stan(model_code = model_string, data = stan_data,iter = 10, chains = 4, seed=seed)    
  }
  
  if(parallel==TRUE){
    # set up for parallel sampling:
    require(doMC)
    options(cores=n_cores)
    registerDoMC()
    
    #prep for parallel case
    fun <- function(i){
      fit <- stan(fit = initial_fit, data = stan_data, iter = iter, warmup = warmup, chains = 1, thin=thin, chain_id=i, seed=seed)
      return(fit)
    }
    
    parallel_fit <- foreach(i = 1:num_chains) %dopar% fun(i)
    
    # squish fit objects together:
    fit <- sflist2stanfit(parallel_fit)
    
  } else {
    fit <- stan(fit = initial_fit, data = stan_data, iter = iter, warmup = warmup, chains = num_chains, thin=thin, seed=seed)
    
  }
  
  return(fit) 
  
}