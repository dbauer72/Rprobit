#' Provide Residuals for a Fitted Model
#'
#' @description
#' This function provides the residuals for a fitted model.
#'
#' @inheritParams AIC.Rprobit_cl
#'
#' @return
#' A \code{matrix}, the residuals.
#'
#' @exportS3Method
#'
#' @importFrom stats residuals

residuals.Rprobit_cl <- function(object, ...) {
  object$residuals()
}
