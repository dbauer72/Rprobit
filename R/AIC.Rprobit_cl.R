#' Calculate the AIC Value Based on a Fitted Model
#'
#' @description
#' This function calculates the AIC value based on a fitted model.
#'
#' @param object
#' An \code{\link{Rprobit_cl}} object.
#' @param ...
#' Not used.
#'
#' @return
#' A \code{numeric}, the AIC value.
#'
#' @exportS3Method 
#'
#' @importFrom stats AIC

AIC.Rprobit_cl <- function(object, ...) {
  object$AIC()
}
