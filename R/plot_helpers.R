# Functions that act as plot helpers for psybayes-type things...

# Plot mean and HDI as ggplot2 pointrange ---------------------------------------------------
#' Mean and HDI of a vector as ggplot2 pointrange.
#' 
#' This is a wrapper for ggplot2's stat_summary function that acts on 
#' a y vector from your ggplot2 object, producing a pointrange geom 
#' that you add to a plot. The point shows the mean of the distribution, 
#' the range shows the 95% HDI.
#' 
#' @export  
#' 
#' @param ... additional values passed to stat_summary.
#' @author Thomas Wallis
#' @examples
#' dat <- data.frame(x = c("Distribution 1", "Distribution 2"),
#' y = rnorm(100))
#' dat$y[dat$x == "Distribution 1"] <- dat$y[dat$x == "Distribution 1"] + 1
#' library(ggplot2)
#' fig <- ggplot(dat, aes(x = x, y = y)) + plot_hdi_pointrange()
#' fig

plot_hdi_pointrange <- function(...){
  require(ggplot2)
  
  hdi_min <- function(x){
    hdis <- hdi(x)
    return(hdis[1])
  }
  
  hdi_max <- function(x){
    hdis <- hdi(x)
    return(hdis[2])
  }
  
  stat_summary(fun.y = mean,
               fun.ymin = hdi_min,
               fun.ymax = hdi_max, 
               geom = "pointrange", 
               ...)
}