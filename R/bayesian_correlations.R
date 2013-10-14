
# Bayesian correlations----------------------------------

# see /Bayesian_Analysis/stan/correlation_basic for testing.

## function to estimate Pearson's R and Spearman's rho between the columns in a data frame.
# input should be a data frame with rows as observations and columns as variables whose correlation structure you wish to estimate.
# all variables must be continuous (no factors).

bayesian_correlations_mcmc <- function(data_frame,parallel=TRUE,iter=1000,n_saved_samples=1000){
  library(rstan)
  
  # Stan data prep -----------------------------------------------------------
  
  scaled <- scale(data_frame)
  fun <- function(i,data_frame){
    return(rank(data_frame[,i]))
  }
  ranked <- sapply(1:ncol(data_frame),fun,data_frame=data_frame)
  ranked <- scale(ranked)
  
  stanDat <- list(scaled = scaled,
                  ranked = ranked,
                  n_obs = nrow(data_frame),
                  n_vars = ncol(data_frame))
  
  # Stan model -----------------------------------------------------------
  model <- '
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
}

transformed parameters {
  cov_matrix[n_vars] Sigma_r;
  cov_matrix[n_vars] Sigma_rho;
  
  # since data are pre-scaled to mean 0, sd 1, dont need to estimate sigma separately... See Stan model p. 95 (Fully Bayes correlated topic model):
  for (m in 1:n_vars){
    Sigma_r[m,m] <- sigma[m] * sigma[m] * Pearson_r[m,m];
    Sigma_rho[m,m] <- sigma[m] * sigma[m] * Spearman_rho[m,m];
    for (n in (m + 1):n_vars){
      Sigma_r[m,n] <- sigma[m] * sigma[n] * Pearson_r[m,n];
      Sigma_r[n,m] <- Sigma_r[m,n];

      Sigma_rho[m,n] <- sigma[m] * sigma[n] * Spearman_rho[m,n];
      Sigma_rho[n,m] <- Sigma_rho[m,n];      
    }
  }
}

model {
  # priors over correlation matrix:
  Pearson_r ~ lkj_corr(1.0);
  Spearman_rho ~ lkj_corr(1.0);

  for (i in 1:n_obs){
    col(scaled_t,i) ~ multi_normal(mu,Sigma_r);
    col(ranked_t,i) ~ multi_normal(mu,Sigma_rho);
  }
}

'
  # do sampling ------------------
  source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')
  
  fit <- stan_sample(model_string = model, stan_data = stanDat, iter = iter, n_saved_samples=n_saved_samples, parallel=parallel)
  
  return(fit)
}


# Plot correlation coefficient histograms for every pairing ---------------------------

bayesian_correlations_plot <- function(data_frame,fit,out_file='correlations.pdf',true_correlations=FALSE){
  library(rstan)
  library(ggplot2)
  
  params <- extract(fit)
  
  # Compute standard correlations -----------------------------------------------------------
  # standard correlation estimates:
  standard_pearson <- cor(data_frame,method='pearson')
  standard_spearman <- cor(data_frame,method='spearman')
  
  # Visually compare correlations with their true values and the standard estimates -----------------------------------------------------------
  
  metric_names <- c("Pearson's r","Spearman's rho")
  
  pdf(out_file)
  
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
      
      fig <- ggplot(markers,aes(x=x,colour=standard)) + facet_wrap(~ metric, ncol=1)
      fig <- fig + geom_histogram(data=sample_frame,aes(x=sample,colour=NULL))
      fig <- fig + geom_vline(data=markers,aes(xintercept=x,colour=standard))
      fig <- fig + xlab('Estimated correlation coefficient') + ylab('Frequency')
      fig <- fig + scale_color_discrete(name='')
      fig <- fig + ggtitle(paste('Correlations for variables ',i,' and ',j,sep=""))
      fig <- fig + theme_grey(base_size=11)
      print(fig)      
    }
  }
  dev.off()
}


# Output a table of centre, HDIs for confidence intervals ---------------------------

bayesian_correlations_table <- function(data_frame,fit,variable_names=FALSE){
  library(rstan)
  library(ggplot2)
  source('~/Dropbox/R_Functions/mcmc_helper_funcs.R')
  
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
        y <- mean(params$Pearson_r[,i,j])
        hdis <- hdiOfMCMC(params$Pearson_r[,i,j])
        
        y <- round(y,digits=2)
        hdis <- round(hdis,digits=2)
        
        # format mean and hdi of this variable into a text string, paste into table:
        text <- paste(y,' [',hdis[1],',',hdis[2],']',sep="")
        pearson_table[i,j] <- text
        
        # for Spearman -----------------------------
        # compute mean and HDI of this cell:
        y <- mean(params$Spearman_rho[,i,j])
        hdis <- hdiOfMCMC(params$Spearman_rho[,i,j])
        
        y <- round(y,digits=2)
        hdis <- round(hdis,digits=2)
        
        # format mean and hdi of this variable into a text string, paste into table:
        text <- paste(y,' [',hdis[1],',',hdis[2],']',sep="")
        spearman_table[i,j] <- text
      }
    }
  }
  
  text <- '1'
  pearson_table[n_vars,n_vars] <- text
  spearman_table[n_vars,n_vars] <- text
  
  return(list(pearson_table=pearson_table,spearman_table=spearman_table))
  
  # this can now be printed using xtable.
}

