#' Adjust Variance Matrix Estimates
#'
#' @description
#' This function adjusts variance matrix estimate if it contains negative
#' eigenvalues.
#'
#' @param mat
#' A variance \code{matrix}.
#' @param tol
#' A \code{numeric}, positive lower bound of acceptable eigenvalues.
#' By default, \code{tol = 1e-6}.
#'
#' @return
#' A \code{list} containing the following elements:
#' * \code{mat}: adjusted \code{matrix}
#' * \code{ok}: a boolean, if \code{FALSE} matrix is adjusted
#' * \code{mat_inv}: inverse of \code{mat}
#'
#' @keywords internal

adjust_variance_matrix <- function(mat, tol = 1e-6) {

  ### make sure that matrix is symmetric
  mat <- 0.5 * (mat + t(mat))

  ### calculate eigenvalues and vectors
  dd <- eigen(mat, symmetric = TRUE)
  lambda <- dd$values

  ### 'ok' flags, if adjustment needs to be made
  ok <- ifelse(min(lambda) < tol, FALSE, TRUE)

  ### adjust eigenvalues and variance matrix
  lambda <- pmax(lambda, tol)
  mat <- dd$vectors %*% diag(lambda) %*% t(dd$vectors)
  mat_inv <- dd$vectors %*% diag(1 / lambda) %*% t(dd$vectors)

  ### return
  list(
    mat = mat,
    ok = ok,
    mat_inv = mat_inv
  )
}
