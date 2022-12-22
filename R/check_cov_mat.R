#' Checks if theta vector creates valid covariance matrices
#'
#' @description
#' This function checks if a valid covariance matrix is created.
#'
#' @param theta
#' A \code{numeric} vector.
#' @param mod
#' A \code{\link{mod_cl}} object.
#'
#' @return
#' Either \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal

check_cov_mat <- function(theta, mod) {
  par <- build_par_from_mod(theta = theta, mod = mod)
  out <- TRUE
  if (mod$lRE > 0) {
    out <- out & (min(eigen(par$Omega[1:mod$lRE, 1:mod$lRE])$values) > 0)
  }
  out <- out & (min(eigen(par$Sigma)$values) >= 0)
  return(out)
}
