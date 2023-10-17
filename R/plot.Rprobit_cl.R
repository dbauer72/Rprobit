#' Plot the coefficient distribution function 
#' (makes sense, if random coefficients or latent classes exist) 
#'
#' @description
#' This function provides a plot for the 
#'
#' @param x 
#' Rprobit_cl object containing the model.
#' @param ...
#' passed on to the plot function of the Rprobit_cl class. 
#' @return
#' A plot.
#'
#' @exportS3Method 

plot.Rprobit_cl <- function(x, ...) {
  x$plot(...)
}
