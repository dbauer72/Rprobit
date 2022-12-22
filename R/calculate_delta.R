#' Build a covariance matrix and calculates the standard deviations
#' @description Function that builds a covariance matrix and calculates the standard deviations using the Delta method.
#' @param HL Description
#' @param fL Description
#' @param thL Description
#' @param alt Number of alternatives
#' @param vv Covariance of estimates
#' @return A list containing the following elements:
#' \item{Sigma}{Built matrix}
#' \item{Sigma_sd}{Standard deviations}
#' @keywords internal

calculate_delta <- function(HL, fL, thL, alt, vv = NULL) {
  n <- dim(HL)[1]
  if (n == 1) {
    Sigma_chol <- HL %*% thL + fL
    Sigma <- matrix(Sigma_chol^2, 1, 1)
    if (!is.null(vv)) {
      var_L <- 4 * Sigma_chol^2 * vv
      sd_Sig <- matrix(sqrt(var_L), 1, 1)
    }
  }
  ### calculate Sigma matrix.
  if (n > 1) {
    eL <- matrixcalc::elimination.matrix(alt)
    Sigma_chol_vech <- HL %*% thL + fL
    Sigma_chol <- matrix(t(eL) %*% Sigma_chol_vech, alt, alt)
    Sigma <- Sigma_chol %*% t(Sigma_chol)
    sd_Sig <- NULL

    if (!is.null(vv)) {
      ### calculate Delta rule
      J_ch_Sig <- J_chol(Sigma_chol_vech)
      vSig_vech <- J_ch_Sig %*% HL %*% vv %*% t(J_ch_Sig %*% HL)
      vSig_vech_sd <- sqrt(diag(vSig_vech))
      D <- duplmat(alt)
      sd_Sig <- matrix(D %*% (vSig_vech_sd), alt, alt)
    }
  }
  ### prepare output with additional sd matrices.
  out <- list(
    "Sigma" = Sigma,
    "Sigma_sd" = sd_Sig
  )
  return(out)
}
