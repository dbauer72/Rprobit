#' Build model parameters from mod-object and compute standard errors
#'
#' @description
#' This function builds model parameters from a \code{mod} object and computes
#' standard errors using the Delta method.
#'
#' @details
#' TODO: The Delta method...
#'
#' @param theta
#' A \code{numeric} vector containing model parameters.
#' @param mod
#' An object of class \code{\link{mod_cl}}.
#' @param variances
#' A \code{matrix} of parameter variances.
#' By default, \code{variances = NA}, in which case no standard errors are
#' calculated.
#' @param labels
#' A \code{list} with two elements:
#' - \code{vars}, a \code{list} of \code{characters}, which are names for the
#'                regressor variables,
#' - \code{cats}, a \code{list} of \code{characters}, which are names for the
#'                categories.
#'
#' By default, \code{labels = NULL}, in which case the regressors and categories
#' are not named.
#'
#' @return A \code{list} containing the following elements:
#' \item{b}{b}
#' \item{b_sd}{b_sd}
#' \item{Omega_chol}{Omega_chol}
#' \item{Omega}{Omega}
#' \item{Omega_sd}{Omega_sd}
#' \item{Sigma_chol}{Sigma_chol}
#' \item{Sigma}{Sigma}
#' \item{Sigma_sd}{Sigma_sd}
#'
#' @keywords internal
#'
#' @importFrom matrixcalc elimination.matrix

build_par_from_mod <- function(theta, mod, variances = NA, labels = NULL) {
  
  ### build 'b' and calculate standard errors
  Hb <- mod$Hb
  fb <- mod$fb
  if (mod$lthb > 0){
    b <- Hb %*% theta[1:mod$lthb] + fb
    if (!any(is.na(variances))) {
      vb <- Hb %*% variances[1:mod$lthb, 1:mod$lthb] %*% t(Hb)
      b_sd <- sqrt(diag(vb))
    } else {
      b_sd <- NULL
    }
  } else {
    b <- fb
    b_sd = NULL
  }
  
  ### build 'Omega' and calculate standard errors
  ### (0 entries in 'Omega' for coefficients without random effect)
  lRE <- mod$lRE
  HO <- mod$HO
  fO <- mod$fO
  Omega_chol <- matrix(0, length(b), length(b))
  if (mod$lthO > 0) {
    Omega_chol_vech <- HO %*% theta[(mod$lthb + 1):(mod$lthb + mod$lthO)] + fO
  } else {
    Omega_chol_vech <- fO
  }
  if (lRE > 1) {
    ### auxiliary matrix to go from vech(Omega) to vec(Omega)
    eO <- matrixcalc::elimination.matrix(lRE)
    Omega_chol[1:lRE, 1:lRE] <- matrix(t(eO) %*% (Omega_chol_vech), lRE, lRE)
  } else if (lRE == 1) {
    Omega_chol[1, 1] <- Omega_chol_vech
  }
  Omega <- Omega_chol %*% t(Omega_chol)
  if (!any(is.na(variances)) & (mod$lthO > 0)) {
    ind_Om_ch <- (mod$lthb + 1):(mod$lthb + mod$lthO)
    sdom <- calculate_delta(HL = mod$HO, fL = mod$fO, thL = theta[ind_Om_ch], alt = lRE, vv = variances[ind_Om_ch, ind_Om_ch])
    Omega_sd <- matrix(0, length(b), length(b))
    Omega_sd[1:lRE, 1:lRE] <- sdom$Sigma_sd
  } else {
    Omega_sd <- NULL
  }
  
  ### build 'Sigma' and calculate standard errors
  HL <- mod$HL
  fL <- mod$fL
  if (mod$ordered) {
    nL <- round(-0.5 + sqrt(2 * length(fL) + 0.25))
  } else {
    nL <- mod$alt
  }
  if (nL>1){
    eL <- matrixcalc::elimination.matrix(nL)
    if (mod$lthL > 0) {
      Sigma_chol_vech <- HL %*% theta[(mod$lthb + mod$lthO + 1):(mod$lthb + mod$lthO + mod$lthL)] + fL
      Sigma_chol <- matrix(t(eL) %*% (Sigma_chol_vech), nL, nL)
    } else {
      Sigma_chol <- matrix(t(eL) %*% (fL), nL, nL)
    }
    Sigma <- Sigma_chol %*% t(Sigma_chol)
    if (!any(is.na(variances)) & (mod$lthL > 0)) {
      ind_Sig_ch <- (mod$lthb + mod$lthO + 1):(mod$lthb + mod$lthO + mod$lthL)
      sdsig <- calculate_delta(HL = mod$HL, fL = mod$fL, thL = theta[ind_Sig_ch], alt = nL, vv = variances[ind_Sig_ch, ind_Sig_ch])
      Sigma_sd <- sdsig$Sigma_sd
    } else {
      Sigma_sd <- NULL
    }
  } else {
    if (mod$lthL > 0) {
      Sigma_chol_vech <- HL %*% theta[(mod$lthb + mod$lthO + 1):(mod$lthb + mod$lthO + mod$lthL)] + fL
      Sigma_chol <- matrix(Sigma_chol_vech, 1, 1)
    } else {
      Sigma_chol <- matrix(fL, 1, 1)
    }
    Sigma <- Sigma_chol %*% Sigma_chol
    if (!any(is.na(variances)) & (mod$lthL > 0)) {
      ind_Sig_ch <- (mod$lthb + mod$lthO + 1):(mod$lthb + mod$lthO + mod$lthL)
      sdsig <- 2*Sigma_chol_vech[1,1]*sqrt(variances[ind_Sig_ch, ind_Sig_ch])
      Sigma_sd <- sdsig
    } else {
      Sigma_sd <- NULL
    }
  }
  
  ##### add names to be variables
  if (!is.null(labels)) {
    if (length(labels[["vars"]]) == length(b)) {
      rownames(b) <- labels[["vars"]]
      if (mod$ordered == FALSE){
        rownames(Sigma) <- labels[["cats"]]
        colnames(Sigma) <- labels[["cats"]]
      }
      rownames(Omega) <- labels[["vars"]]
      colnames(Omega) <- labels[["vars"]]
      
      if (!is.null(b_sd)) {
        rownames(b_sd) <- labels[["vars"]]
        if (mod$ordered == FALSE){
          rownames(Sigma_sd) <- labels[["cats"]]
          colnames(Sigma_sd) <- labels[["cats"]]
        }
        rownames(Omega_sd) <- labels[["vars"]]
        colnames(Omega_sd) <- labels[["vars"]]
      }
    }
  }
  
  out <- list(
    "b" = b,
    "b_sd" = b_sd,
    "Omega_chol" = Omega_chol,
    "Omega" = Omega,
    "Omega_sd" = Omega_sd,
    "Sigma_chol" = Sigma_chol,
    "Sigma" = Sigma,
    "Sigma_sd" = Sigma_sd
  )
  
  ###### add tauk for ordered case
  if (mod$ordered) {
    theta_tauk <- theta[(mod$lthb + mod$lthO + mod$lthL + 1):(mod$lthb + mod$lthO + mod$lthL + mod$alt-1)]
    if (length(theta_tauk) != mod$alt - 1) {
      message("Length of parameters for tauk does not match mod$alt.")
    } else {
      tauk <- rep(0, mod$alt - 1)
      tauk[1] <- theta_tauk[1]
      dtauk <- c(tauk[1], exp(theta_tauk[-1]))
      tauk <- cumsum(dtauk)
      out$tauk <- tauk
    }
    
    # state space system, if 
  }
  
  return(out)
}
