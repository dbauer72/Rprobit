#' Provide the Log-Likelihood Value for a Fitted Model
#'
#' @description
#' This function provides the log-likelihood value for a fitted model.
#'
#' @inheritParams AIC.Rprobit_cl
#'
#' @return
#' A \code{numeric}, the log-likelihood value.
#'
#' @exportS3Method
#'
#' @importFrom stats logLik

logLik.Rprobit_cl <- function(object, ...) {
  object$ll
}
