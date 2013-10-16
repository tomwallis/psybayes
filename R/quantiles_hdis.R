## functions for examining statistics from vectors and matrices (usually in the context of mcmc output).

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
#' @family distribution_summaries
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
  names(hdi_lim) <- c('hdi_min','hdi_max')
  return( hdi_lim )
}


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
#' @family distribution_summaries
#' @examples
#' X <- matrix(rnorm(1000),ncol=10)
#' col_quantile(X)

col_quantile <- function(X,probs=c(0.025,0.975)){
  n_cols <- NCOL(X)
  fun <- function(i,X,probs) quantile(X[,i],probs=probs)
  q <- sapply(1:n_cols,fun,X,probs=probs)
  return(q)
}



# column-wise HDI ---------------------------
#' Calculate the HDI along the columns of a matrix.
#' 
#' Receives a matrix X as input, and calculates the HDI along the columns using sapply.
#' 
#' @export  
#' 
#' @param X a matrix of numbers
#' @param probs vector of quantiles to compute.
#' @return a matrix where columns correspond to columns of \code{X}, and rows to the lower and upper HDI limits.
#' @author Thomas Wallis
#' @family distribution_summaries
#' @examples
#' X <- matrix(rnorm(1000),ncol=10)
#' col_hdi(X)

col_hdi <- function(X,cred_mass=0.95){
  n_cols <- NCOL(X)
  fun <- function(i,X,cred_mass) hdi(X[,i],cred_mass = cred_mass)
  q <- sapply(1:n_cols, fun, X, cred_mass)
  return(q)
}
