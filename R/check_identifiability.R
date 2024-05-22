#' Test for model identifiability
#'
#' @description
#' This function tests whether a given model specification seems identifiable.
#'
#' @param mod
#' An object of class \code{\link{mod_cl}}.
#' @param tol
#' A positive \code{numeric} to define the check tolerance.
#' By default, \code{tol = 1e-6}.
#'
#' @return
#' Either
#' - \code{TRUE} (model seems identifiable) or
#' - \code{FALSE} (a reason for non-identifiability is found).
#'
#' @export
#'
#' @examples
#' ### this model is not identified (beta_3 = beta_1 + beta_2, all estimated)
#' mod1 <- mod_cl$new(
#'   Hb  = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 0), nrow = 3, ncol = 3),
#'   fb  = matrix(c(0, 0, 0), ncol = 1),
#'   HO  = matrix(0, 0, 0),
#'   fO  = matrix(0, 0, 0),
#'   HL  = matrix(0, nrow = 3, ncol = 0),
#'   fL  = matrix(c(0, 0, 1), ncol = 1),
#'   alt = 2
#' )
#' check_identifiability(mod1)
#'
#' ### not estimating beta_3 makes the model identified
#' mod2 <- mod1$clone()
#' mod2$Hb <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, ncol = 2)
#' check_identifiability(mod2)
#'
#' ### this model is not identified (no scale restriction)
#' mod3 <- mod_cl$new(
#'   Hb  = diag(2),
#'   fb  = matrix(c(0, 0), ncol = 1),
#'   HO  = matrix(0, 0, 0),
#'   fO  = matrix(0, 0, 0),
#'   HL  = diag(3),
#'   fL  = as.matrix(c(0, 0, 0), ncol = 1),
#'   alt = 2
#' )
#' check_identifiability(mod3)
#'
#' ### fixing one beta entry ensures identifiability
#' mod4 <- mod3$clone()
#' mod4$Hb <- matrix(c(0, 1), ncol = 1)
#' mod4$fb <- matrix(c(1, 0), ncol = 1)
#' check_identifiability(mod4)
#'
#' ### also fixing one Sigma entry ensures identifiability
#' mod5 <- mod3$clone()
#' mod5$HL <- matrix(0, nrow = 3, ncol = 0)
#' mod5$fL <- matrix(c(0, 0, 1), ncol = 1)
#' check_identifiability(mod5)
#'
#' @importFrom stats rnorm runif

check_identifiability <- function(mod, tol = 1e-6) {
  npar <- mod$lthb + mod$lthL + mod$lthO
  theta <- stats::rnorm(npar)
  if (mod$ordered) {
    theta <- c(theta, -2, stats::runif(mod$alt - 2, -1, 1))
  }
  J <- mod$alt
  P <- dim(mod$Hb)[1]
  dims <- P + P^2 + ifelse(mod$ordered, 1 + J - 2, J^2)
  systems <- matrix(0, dims, npar)
  prepare_vech <- function(theta) {
    par <- build_par_from_mod(theta, mod)
    Sigma <- par$Sigma
    b <- par$b
    norm_b <- norm(b)
    Omega <- par$Omega
    if (!mod$ordered) {
      delta_1 <- diag(J)
      delta_1[, 1] <- -1
      Sigma <- delta_1 %*% Sigma %*% t(delta_1)
      tauk <- numeric()
    } else {
      Sigma <- Sigma[1, 1]
      tauk <- par$tauk[-1] / norm_b
    }
    b <- b / norm_b
    Omega <- Omega / norm_b^2
    Sigma <- Sigma / norm_b^2
    c(as.vector(b), as.vector(Omega), as.vector(Sigma), as.vector(tauk))
  }
  vech_0 <- prepare_vech(theta)
  for (i in 1:npar) {
    theta_shift <- theta
    theta_shift[i] <- theta_shift[i] + 1
    vech_j <- prepare_vech(theta_shift)
    systems[, i] <- (vech_j - vech_0) / tol
  }
  if(abs(det(t(systems) %*% systems)) < tol) {
    warning(
      "Parameter set seems not identifed.\n",
      "Check model specification.",
      call. = FALSE
    )
    return(FALSE)
  }
  theta_scaled <- theta * 2
  vech_scaled <- prepare_vech(theta_scaled)
  if (norm(as.matrix(vech_0 - vech_scaled)) < tol) {
    warning(
      "Scale does not change system.\n",
      "Consider fixing one entry of beta or Sigma.",
      call. = FALSE
    )
    return(FALSE)
  }
  return(TRUE)
}
