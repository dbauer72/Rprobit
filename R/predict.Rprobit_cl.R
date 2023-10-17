#' Prediction Based on a Fitted Model
#' 
#' @description 
#' This function calculates the in-sample predictions for a fitted model.
#' 
#' @inheritParams AIC.Rprobit_cl
#' @param newdata 
#' Either \code{\link{data_raw_cl}} object. Or a data.frame of the same structure as the ones on which 
#' the \code{\link{Rprobit_cl}} \code{object} is based upon.  
#' 
#' @return 
#' A \code{list} of matrices with information on the estimated model.
#' 
#' @exportS3Method 

predict.Rprobit_cl <- function(object, newdata = NULL, ...) {
  object$predict(newdata,...)
}
