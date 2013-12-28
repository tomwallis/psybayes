
# Do mcmc sampling ----------------------------------
#' Estimate Bayesian (posterior) correlation coefficients. 
#' 
#' Estimate Pearson's r and Spearman's rho between the columns in a data frame.
#' This is a wrapper function for Stan.
#' 
#' @export  
#' @import rstan
#' 
#' @param data_frame    A data frame where rows are observations and columns 
#' are the variables whose correlation structure you wish to estimate. All variables must be continuous (no factors).
#' @param prior         A text string containing the (Stan model) prior specification for the correlation matrices.
#' @param ...           Additional variables passed to stan_sample function.
#' 
#' @return Returns an rstan fit object.
#' 
#' @details
#' This function takes a numerical data matrix X as input. From this, the Pearson r and Spearman's rho
#' (rank order correlation coefficient) will be estimated using MCMC sampling from the \link{stan} package.
#' Currently the function standardises the data internally (using the \link{scale} function), but could be easily
#' modified to also estimate the mean and variance of each variable.
#' 
#' The \code{prior} is applied to each correlation matrix (Pearson and Spearman). The default is a multivariate 
#' version of a Beta distribution with parameters (1,1), i.e. it is multivariate uniform, scaled over the range 
#' -1 to 1. See Stan manual for details.
#' 
#' For a big example of all this in action, run \code{example(bcor_mcmc)}.
#' 
#' Note: on my todo list is to fix the warning output about incrementing the log Jacobian. 
#' However, this doesn't seem to be a problem for the sampling. See \url{https://groups.google.com/forum/#!topic/stan-users/9GCAlBWVBtI}.
#' "Like some of our other warning messages, there are "false positives" 
#' from this because some of the transforms applied 
#' (like taking the column of a matrix using function col()) have 
#' a zero log Jacobian adjustment. "
#' 
#' 
#' @author Thomas Wallis.
#' @seealso \link{stan_sample}, \link{stan}
#' @family bcor
#' @examples
#' # generate random variables with a specified correlation structure (from: http://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/)
#' # n_var variables. Assume that they are standardised (mean 0 sd 1).
#' n_vars <- 4
#' # number of observations to simulate
#' n_obs = 100
#' # correlation matrix: 
#' M = matrix(c(1, 0.7, 0.7, 0.5,
#'              0.7, 1, 0.95, 0.3,
#'              0.7, 0.95, 1, 0.3,
#'              0.5, 0.3, 0.3, 1), nrow=4, ncol=4)
#'              
#' # Cholesky decomposition
#' L <- chol(M)
#' nvars = dim(L)[1]            
#' # Random variables that follow an M correlation matrix
#' r <- t(L) %*% matrix(rnorm(n_vars*n_obs), nrow=n_vars, ncol=n_obs)
#' r <- t(r)
#' dat <- data.frame(r)
#' 
#' # Now that we've generated data, fit and analyse:
#' fit <- bcor_mcmc(dat)
#' print(fit)
#' bcor_plot(dat,fit,true_correlation=M)
#' tables <- bcor_table(dat,fit)
#' 
#' # try with a weakly-informative prior (more density over unit diagonal -- should be closer to zero with less data)
#' fit_2 <- bcor_mcmc(dat, prior = 'lkj_corr(2.0)')
#' print(fit_2)
#' bcor_plot(dat,fit_2,true_correlation=M)
#' tables_2 <- bcor_table(dat,fit_2)
#' 
#' tables$pearson_table
#' tables_2$pearson_table
#' cor(dat)

bcor_mcmc <- function(data_frame, prior = 'lkj_corr(1.0)', ...){
  require(rstan)
  
  # Stan data prep -----------------------------------------------------------
  
  scaled <- scale(data_frame)
  fun <- function(i,data_frame){
    return(rank(data_frame[,i]))
  }
  ranked <- sapply(1:ncol(data_frame),fun,data_frame=data_frame)
  ranked <- scale(ranked)
  
  stan_data <- list(scaled = scaled,
                  ranked = ranked,
                  n_obs = nrow(data_frame),
                  n_vars = ncol(data_frame))
  
  # Stan model -----------------------------------------------------------
  model <- paste0('
  data{
  int<lower=1> n_obs;
  int<lower=2> n_vars;
  matrix[n_obs,n_vars] scaled;
  matrix[n_obs,n_vars] ranked;
  }
  
  transformed data {
  # these could easily be estimated from the data, rather than pre-scaling. Just not sure how to rank in Stan (to calculate Spearman).
  vector[n_vars] mu;    
  vector[n_vars] sigma; 
  matrix[n_vars,n_obs] scaled_t;
  matrix[n_vars,n_obs] ranked_t;
  
  for (i in 1:n_vars){
  mu[i] <- 0.0;
  sigma[i] <- 1.0;
  }
  
  # transpose data matrices for indexing with multi_norm, which seems to want column vectors:
  scaled_t <- scaled\';
  ranked_t <- ranked\';
}

parameters {
  corr_matrix[n_vars] Pearson_r;
  corr_matrix[n_vars] Spearman_rho;
  corr_matrix[n_vars] Prior_samples;
}

transformed parameters {
  cov_matrix[n_vars] Sigma_r;
  cov_matrix[n_vars] Sigma_rho;
  cov_matrix[n_vars] Sigma_prior;
  
  # since data are pre-scaled to mean 0, sd 1, dont need to estimate sigma separately... See Stan model p. 95 (Fully Bayes correlated topic model):
  for (m in 1:n_vars){
    Sigma_r[m,m] <- sigma[m] * sigma[m] * Pearson_r[m,m];
    Sigma_rho[m,m] <- sigma[m] * sigma[m] * Spearman_rho[m,m];
    Sigma_prior[m,m] <- sigma[m] * sigma[m] * Prior_samples[m,m];

    for (n in (m + 1):n_vars){
      Sigma_r[m,n] <- sigma[m] * sigma[n] * Pearson_r[m,n];
      Sigma_r[n,m] <- Sigma_r[m,n];

      Sigma_rho[m,n] <- sigma[m] * sigma[n] * Spearman_rho[m,n];
      Sigma_rho[n,m] <- Sigma_rho[m,n];   

      Sigma_prior[m,n] <- sigma[m] * sigma[n] * Prior_samples[m,n];
      Sigma_prior[n,m] <- Sigma_prior[m,n];
    }
  }
}

model {
  # priors over correlation matrix:
  Pearson_r ~ ', prior,';
  Spearman_rho ~ ', prior, ';
  Prior_samples ~ ', prior,';

  for (i in 1:n_obs){
    col(scaled_t,i) ~ multi_normal(mu,Sigma_r);
    col(ranked_t,i) ~ multi_normal(mu,Sigma_rho);
    # prior matrix is not updated with data.
  }
}

')
  
  fit <- stan_sample(model_string = model, data = stan_data, ...)
  
  return(fit)
}


# Produce plots ---------------------------
#' Plot posterior correlation coefficients. 
#' 
#' Plot the correlation coefficients estimated by bcor_mcmc.
#' 
#' @export  
#' @import ggplot2
#' 
#' @param data_frame    A data frame where rows are observations and columns 
#' are the variables whose correlation structure you wish to estimate. 
#' All variables must be continuous (no factors).
#' @param fit           The stanfit object produced by \link{bcor_mcmc}.
#' @param out_file      The filename of the output file (with .pdf extension), or 'none' (print to R directly).
#' @param true_correlations   If a matrix of numbers is passed here, 
#' this is taken to be the True correlation coefficients used to generate the data. 
#' If left at FALSE, the generating correlation coefficients are not plotted.
#' @param plot_prior    If TRUE, plot an estimate of the prior density from the model.
#' 
#' @return If \code{out_file} contains a filename with .pdf extension, 
#' a plot of each pairwise correlation is sent to out_file.
#' Otherwise these are printed to R. 
#' The dark histogram shows samples from the posterior, the faint histogram
#' shows samples from the prior.
#' @author Thomas Wallis.
#' @family bcor
#' @examples
#' See the example for bcor_mcmc.

bcor_plot <- function(data_frame,fit, out_file='none', true_correlations=FALSE, plot_prior=TRUE){
  require(rstan)
  require(ggplot2)
  
  params <- extract(fit)
  
  # Compute standard correlations -----------------------------------------------------------
  # standard correlation estimates:
  standard_pearson <- cor(data_frame,method='pearson')
  standard_spearman <- cor(data_frame,method='spearman')
  
  # Visually compare correlations with their true values and the standard estimates -----------------------------------------------------------
  
  metric_names <- c("Pearson's r","Spearman's rho")
  
  if(out_file!='none') pdf(out_file)
  
  for (i in 1:(dim(data_frame)[2]-1)){
    for (j in (i+1):dim(data_frame)[2]){
      sample_frame <- expand.grid(sample=rep(NA,length=length(params$lp__)),
                                  metric=metric_names)
      
      sample_frame$sample[sample_frame$metric==metric_names[1]] <- params$Pearson_r[,i,j]
      sample_frame$sample[sample_frame$metric==metric_names[2]] <- params$Spearman_rho[,i,j]
      
      if(is.logical(true_correlations)){
        markers <- expand.grid(x=NA,metric=metric_names,standard=c("Standard estimate","Posterior Mean"))
      } else {
        markers <- expand.grid(x=NA,metric=metric_names,standard=c("True","Standard estimate","Posterior Mean"))
        markers$x[markers$metric==metric_names[1] & markers$standard=="True"] <- true_correlations[i,j]
        markers$x[markers$metric==metric_names[2] & markers$standard=="True"] <- true_correlations[i,j]        
      }
      markers$x[markers$metric==metric_names[1] & markers$standard=="Standard estimate"] <- standard_pearson[i,j]
      markers$x[markers$metric==metric_names[2] & markers$standard=="Standard estimate"] <- standard_spearman[i,j]
      markers$x[markers$metric==metric_names[1] & markers$standard=="Posterior Mean"] <- mean(params$Pearson_r[,i,j])
      markers$x[markers$metric==metric_names[2] & markers$standard=="Posterior Mean"] <- mean(params$Spearman_rho[,i,j])
      
      fig <- ggplot(markers,aes(x=x,colour=standard)) + facet_wrap(~ metric, ncol=2)
      
      if(plot_prior==TRUE){
        prior_frame <- data.frame(sample=rep(NA,length=length(params$lp__)))
        prior_frame$sample <- params$Prior_samples[,i,j]
        fig <- fig + stat_bin(data=prior_frame,aes(x=sample,y=..density.., colour=NULL), fill = "black", alpha = 0.25, binwidth=0.03)
      }      
      
      fig <- fig + stat_bin(data=sample_frame,aes(x=sample, y = ..density.., colour=NULL), binwidth=0.03)
      fig <- fig + geom_vline(data=markers,aes(xintercept=x,colour=standard))
      fig <- fig + scale_color_discrete(name='')
      fig <- fig + scale_x_continuous(limits=c(-1,1))
      fig <- fig + xlab('Estimated correlation coefficient') + ylab('Normalised density')
      fig <- fig + ggtitle(paste('Correlations for variables ',i,' and ',j,sep=""))
      fig <- fig + theme_grey(base_size=11)
      print(fig)      
    }
  }
  if(out_file!='none') dev.off()
}


# Produce tables ---------------------------
#' Table of correlation coefficients with HDIs.
#' 
#' Generate a numerical matrix of the correlation coefficients estimated by bcor_mcmc.
#' 
#' @export  
#' @import ggplot2
#' 
#' @param data_frame    A data frame where rows are observations and columns 
#' are the variables whose correlation structure you wish to estimate. 
#' All variables must be continuous (no factors).
#' @param fit           The stanfit object produced by \link{bcor_mcmc}.
#' @param variable_names      Either \code{FALSE} (columns in table will be named after data_frame) or
#' a character vector of length = ncol(data_frame) containing a label for each variable.
#' @param central_tendency_function   A function to calculate the data's central tendency. Should return a
#' scalar with input as a vector. 
#' @param spread_function             A function to calculate the data's spread (e.g. hdi or quantiles).
#' Should return two values (lower, upper) given input as a vector.
#' 
#' @return A list of tables suitable for passing to the xtable function. One table is for the Pearson 
#' Product Moment correlation coefficent, one is for Spearman's rho (rank order correlation) and the 
#' third are the same computed from samples of the prior distribution.
#' @author Thomas Wallis.
#' @family bcor
#' @examples
#' See the example for bcor_mcmc.

bcor_table <- function(data_frame,fit,variable_names=FALSE,
                       central_tendency_function = mean,
                       spread_function = hdi){
  require(rstan)
  
  params <- extract(fit)
  n_samples <- length(params$lp__)
  n_vars <- dim(params$Pearson_r)[2]
  
  pearson_table <- matrix(rep('--',times=n_vars^2),ncol=n_vars)
  
  if (is.logical(variable_names)){
    colnames(pearson_table) <- names(data_frame)
    rownames(pearson_table) <- names(data_frame)
  } else {
    colnames(pearson_table) <- variable_names
    rownames(pearson_table) <- variable_names
  }
  
  spearman_table <- pearson_table
  prior_table <- pearson_table
  
  # this is quite shitty, but I think the best way to look at this is by just printing each cell into a table then using xtable...
  for(i in 1:(n_vars-1)){
    for (j in i:n_vars){
      if (i == j){
        text <- '1'
        pearson_table[i,j] <- text
        spearman_table[i,j] <- text
      } else {
        # for Pearson -----------------------------
        # compute mean and HDI of this cell:
        y <- central_tendency_function(params$Pearson_r[,i,j])
        hdis <- spread_function(params$Pearson_r[,i,j])
        
        y <- round(y,digits=2)
        hdis <- round(hdis,digits=2)
        
        # format mean and hdi of this variable into a text string, paste into table:
        text <- paste(y,' [',hdis[1],',',hdis[2],']',sep="")
        pearson_table[i,j] <- text
        
        # for Spearman -----------------------------
        # compute mean and HDI of this cell:
        y <- central_tendency_function(params$Spearman_rho[,i,j])
        hdis <- spread_function(params$Spearman_rho[,i,j])
        
        y <- round(y,digits=2)
        hdis <- round(hdis,digits=2)
        
        # format mean and hdi of this variable into a text string, paste into table:
        text <- paste(y,' [',hdis[1],',',hdis[2],']',sep="")
        spearman_table[i,j] <- text
        
        # for Priors -----------------------------
        # compute mean and HDI of this cell:
        y <- central_tendency_function(params$Prior_samples[,i,j])
        hdis <- spread_function(params$Prior_samples[,i,j])
        
        y <- round(y,digits=2)
        hdis <- round(hdis,digits=2)
        
        # format mean and hdi of this variable into a text string, paste into table:
        text <- paste(y,' [',hdis[1],',',hdis[2],']',sep="")
        prior_table[i,j] <- text
      }
    }
  }
  
  text <- '1'
  pearson_table[n_vars,n_vars] <- text
  spearman_table[n_vars,n_vars] <- text
  prior_table[n_vars,n_vars] <- text
  
  return(list(pearson_table=pearson_table,spearman_table=spearman_table,prior_table=prior_table))
  
  # this can now be printed using xtable.
}

