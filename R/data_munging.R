## functions for arranging data.
# Binomial design matrix ----------------------------------
#' Binomial design matrix from Bernoulli trials.
#' 
#' Transform a data frame where each row is a Bernoulli trial into 
#' a binomial data frame (ntrials, nsucceses) based on the design matrix for a model.
#' 
#' @export  
#' 
#' @import plyr
#' 
#' @param data     A data frame with rows as Bernoulli trials.
#' @param formula  The formula for the model design at the lowest level of the model.
#' @param group_factor  A character vector specifying the name of one or more grouping variables found in \code{data}. Usually "subject". 
#' @param response  A character vector specifying the name of the response variable (correct or incorrect).
#' 
#' @details
#' In Stan it is often much more efficient to sample a model binomially 
#' rather than as individual Bernoulli trials. This function can be used
#' to transform a data frame of Bernoulli trials (like you might save in an 
#' experiment) into a binomial data frame, based on a formula for the model.
#' 
#' @return Returns a new data frame.
#' 
#' @family data_munging
#' 
#' @author Thomas Wallis.
#' @examples
#' data(detection)
#' head(detection)
#' # detection has rows as Bernoulli trials.
#' # we would like to model performance as a function of log contrast, for each subject
#' # within a population.
#' binom_list <- bern2binom(detection, correct ~ log(contrast), group_factors = "subject", response = "correct")
#' 
#' binom_list$binom_dat
#' 
#' # design matrix can be passed to Stan along with response and grouping variables:
#' binom_list$design_matrix

bern2binom <- function(data, formula, group_factors='none', response='y'){
  require(plyr)
  
  # prepare model matrix:
  x <- model.matrix(as.formula(formula), data=data)
  x <- data.frame(x)
  beta_names <- colnames(x)
  
  # add response variable:
  x$y <- data[,response]
  
  # add grouping variable if specified:
  if(group_factors != 'none' | length(group_factors) > 1){
    for(i in 1:length(group_factors)){
      text <- paste0('x$',group_factors[i],' <- data$',group_factors[i])
      eval(parse(text=text))
    }
    # create binomial:
    binom_dat <- ddply(x,(c(beta_names,group_factors)),summarise,n_cor=sum(y),n_trials=length(y))
  } else {
    # create binomial:
    binom_dat <- ddply(x,(beta_names),summarise,n_cor=sum(y),n_trials=length(y))
  }
  
  
  # create design matrix again by removing other variables:
  X <- binom_dat[,beta_names]
  
  # return a list of both the new matrices:
  return(list(binom_dat = binom_dat, design_matrix = X))
}
