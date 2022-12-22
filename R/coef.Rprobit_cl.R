#' Coefficients of a Fitted Probit Model
#' 
#' @description 
#' This function returns the regression coefficients for a fitted probit model.
#' 
#' @inheritParams AIC.Rprobit_cl
#' 
#' @return 
#' A \code{numeric} vector of coefficients.
#' 
#' @exportS3Method 
#' 
#' @importFrom stats coef

coef.Rprobit_cl <- function(object, ...) {
  object$coef()
}
