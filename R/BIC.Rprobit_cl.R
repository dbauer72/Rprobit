#' Calculate the BIC Value Based on a Fitted Model
#'
#' @description
#' This function calculates the BIC value based on a fitted model.
#'
#' @inheritParams AIC.Rprobit_cl
#'
#' @return
#' A \code{numeric}, the BIC value.
#'
#' @exportS3Method 
#'
#' @importFrom stats BIC

BIC.Rprobit_cl <- function(object, ...) {
  object$BIC()
}
