#' Prediction Based on a Fitted Model
#' 
#' @description 
#' This function calculates the in-sample predictions for a fitted model.
#' 
#' @inheritParams AIC.Rprobit_cl
#' 
#' @return 
#' A \code{list} of matrices with information on the estimated model.
#' 
#' @exportS3Method 

predict.Rprobit_cl <- function(object, ...) {
  object$predict()
}
