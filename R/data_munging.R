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


# Bin bernoulli trials along a continuous x -------------------------
#' Bin bernoulli trials on a continuous variable x, calculate stats.
#' 
#' Takes an x and a y vector, where y is binary (bernoulli trials).
#' 
#' @export  
#' @import plyr
#' 
#' @return a data frame with x bins as rows, and the number of successes, fails, and trials in each bin.
#' 
#' @param data   a data frame containing (at least) a numeric vector \code{x} and a binary outcome vector \code{y}.
#' @param breaks  a single number giving the number of bins desired, or a vector of breakpoints between cells.
#' @param spacing  a character vector of either \code{"equal"} or \code{"quantile"}, producing equally spaced
#' bins along x or bins with approximately equal numbers of data points, respectively. 
#' Doesn't do anything if you pass \code{breaks} as a vector of breakpoints.
#' @param x_name A character vector specifying the name of the continuous variable to bin over.
#' @param y_name A character vector specifying the name of the binary outcome variable.
#' @param additional_factors  A character vector specifying any factors in the dataset to split x over.
#' @param ... additional arguments passed to \link{beta_cis}.
#'
#' @author Thomas Wallis
#' @seealso \link{beta_cis}
#' @examples
#' library(ggplot2)
#' dat <- data.frame(x = rnorm(100), y = rbinom(100,1,prob=0.5))
#' qplot(dat$x,dat$y)
#' equal_spacing <- bern_bin(dat,breaks=5)
#' equal_numbers <- bern_bin(dat,breaks=5,spacing='quantile') 
#' custom_spacing <- bern_bin(dat,breaks=seq(-3,3,l=6))
#' 
#' # plot:
#' equal_spacing$label <- 'Equal spacing'
#' equal_numbers$label <- 'Equal numbers'
#' custom_spacing$label <- 'Custom'
#' 
#' plot_dat <- rbind(equal_spacing,equal_numbers,custom_spacing)
#' 
#' fig <- ggplot(plot_dat,aes(x=xmid,xmax=xmax,xmin=xmin,y=y,ymax=ymax,ymin=ymin)) + geom_point() + geom_errorbar() + geom_errorbarh()
#' fig <- fig + facet_wrap(~ label, ncol=1)
#' fig
#' 
#' ## Show use of additional factors, and changing beta distribution calculations:
#' dat$a_factor <- rbinom(nrow(dat),1,prob=0.5)
#' dat$a_factor <- factor(dat$a_factor)
#' 
#' binned <- bern_bin(dat,breaks=5, additional_factors = "a_factor", probs = c(0.025, 0.975))
#' fig <- ggplot(binned,aes(x=xmid,xmax=xmax,xmin=xmin,y=y,ymax=ymax,ymin=ymin)) + geom_point() + geom_errorbar() + geom_errorbarh()
#' fig <- fig + facet_wrap(~ a_factor, ncol=1)
#' fig

bern_bin <- function(data, breaks, spacing = 'equal', 
                     x_name = "x",
                     y_name = "y",
                     additional_factors = "none", ...){
  
  # set up key variables:
  text <- paste0('x <- data$',x_name)
  eval(parse(text=text))
  text <- paste0('y <- data$',y_name)
  eval(parse(text=text))
  
  # from http://r.789695.n4.nabble.com/Obtaining-midpoints-of-class-intervals-produced-by-cut-and-table-td908058.html
  cut2num <- function(f){ 
    labs <- levels(f) 
    d <- data.frame(xmin = c(NA, as.numeric( sub("\\((.+),.*", "\\1", labs[-1]) )), 
                    xmax = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
    # correct for lower=TRUE in regexp (causes regexp not to match first label:
    d$xmin[1] <- as.numeric(sub("\\[(.+),.*", "\\1", labs)[1])
    d$xmid <- rowMeans(d) 
    return(d)
#     rm(d)
  } 
  
  produce_df <- function(breaks){
    # use cut to produce factors. Find midpoints using the cut2num function.
    x_cut <- cut(x,breaks=breaks,include.lowest=TRUE)
    cut_limits <- cut2num(x_cut)
    
    if(additional_factors == "none"){
      d <- data.frame(x = x_cut, y = y)
      binomial_df <- ddply(d,.(x),summarise,n_success = sum(y), n_trials = length(y))  
    } else {
      d <- data.frame(x = x_cut, y = y)
      d <- cbind(d, subset(data, select = additional_factors))
      factors <- c('x',additional_factors)
      binomial_df <- ddply(d,factors,summarise,n_success = sum(y), n_trials = length(y))  
    }
    
    instances <- nrow(binomial_df) / nrow(cut_limits)
    df <- do.call("rbind", replicate(instances, cut_limits, simplify = FALSE))
    df <- cbind(df, subset(binomial_df,select = -x))
    df$n_fails <- binomial_df$n_trials - binomial_df$n_success
    df$prop_corr <- binomial_df$n_success / binomial_df$n_trials
    
    df <- beta_cis(df, ...)
    
    return(df)
  }
  
  
  

  if(length(breaks)>1){
    df <- produce_df(breaks)
    return(df)
  } else {
    if(spacing=="equal"){
      df <- produce_df(breaks)
      return(df)
    } 
    
    if(spacing=="quantile") {
      # use quantiles to produce break points.
      breaks <- quantile(x, probs = seq(0, 1, length = breaks + 1))
      df <- produce_df(breaks)
      return(df)
    }
    
    if(spacing != "equal" & spacing != "quantile") stop("improperly specified spacing")
  }
}

# Beta distribution CIs -------------------------
#' From binomial data frame, compute confidence limits of beta distribution.
#' 
#' Input is a data frame containing columns n_success, n_fails.
#' 
#' @export  
#' 
#' @return the data frame with added y, ymin and ymax columns. Data (n_successes, etc) 
#' will be original, not rule-of-succession corrected.
#' 
#' @param data   a data frame.
#' @param probs   probabilities to return (increasing). If this is a length 3 vector, 
#' the values will be assigned to ymin, y and ymax. If it's a vector of 
#' length 2, then this specifies ymin and ymax, and y will be the proportion correct.
#' @param rule_of_succession  If true, add one success and one failure to every observation
#' (Laplace's rule of succession correction).
#'
#'@seealso \link{bern_bin}
#'
#' @author Thomas Wallis
#' @examples
#' example(bern_bin)

beta_cis <- function(data, probs = c(0.025, 0.5, 0.975), rule_of_succession = TRUE){
  d <- data
  
  if(rule_of_succession == TRUE){
    d$n_success <- data$n_success + 1
    d$n_fails <- data$n_fails + 1
    d$n_trials <- data$n_trials + 2
  }
  

  if(length(probs) == 3){
    data$ymin <- qbeta(probs[1],d$n_success,d$n_fails)
    data$ymax <- qbeta(probs[3],d$n_success,d$n_fails)
    data$ymid <- qbeta(probs[2],d$n_success,d$n_fails)
  }
  
  if(length(probs) == 2){
    data$ymin <- qbeta(probs[1],d$n_success,d$n_fails)
    data$ymax <- qbeta(probs[2],d$n_success,d$n_fails)
    data$ymid <- data$n_success / (data$n_success + data$n_fails)
  }
  
  if(length(probs) != 3 & length(probs) != 2) 
    stop("probabilities must be a vector of length 2 or 3")
  
  return(data)
}
