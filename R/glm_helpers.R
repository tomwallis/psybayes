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

