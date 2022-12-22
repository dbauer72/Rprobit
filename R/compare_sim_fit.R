#' Compare true with fitted model
#' @description Function that compares true model parameters of the dgp for a simulated data set with parameters of a fitted model.
#' @param theta_0 True model parameters for simulated data
#' @param Rprobit_obj \code{Rprobit_obj}-object
#' @param show_theta boolean, if TRUE (default) compare theta vector, otherwise compare entries of b, Omega and Sigma
#' @return matrix with theta vectors in rows
#' @keywords internal

compare_sim_fit <- function(theta_0, Rprobit_obj, show_theta = TRUE) {
  if (show_theta) {
    out <- rbind(true = theta_0, estimate = Rprobit_obj$theta)
    colnames(out) <- c(paste(rep("b", Rprobit_obj$mod$lthb), seq_len(Rprobit_obj$mod$lthb), sep = "_"), paste(rep("o", Rprobit_obj$mod$lthO), seq_len(Rprobit_obj$mod$lthO), sep = "_"), paste(rep("l", Rprobit_obj$mod$lthL), seq_len(Rprobit_obj$mod$lthL), sep = "_"))
  } else {
    out <- numeric(0)

    if (Rprobit_obj$mod$lthb > 0) {
      out_b <- rbind(true = c(Rprobit_obj$mod$Hb %*% theta_0[1:Rprobit_obj$mod$lthb] + Rprobit_obj$mod$fb), estimate = c(Rprobit_obj$mod$Hb %*% Rprobit_obj$theta[1:Rprobit_obj$mod$lthb] + Rprobit_obj$mod$fb))
      colnames(out_b) <- paste("b", 1:length(Rprobit_obj$mod$fb), sep = "_")
      out <- cbind(out, out_b)
    }

    if (Rprobit_obj$mod$lthO > 0) {
      if (Rprobit_obj$mod$lRE > 1) {
        Omega_vec_est <- t(matrixcalc::elimination.matrix(Rprobit_obj$mod$lRE)) %*% (Rprobit_obj$mod$HO %*% Rprobit_obj$theta[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fO)
        Omega_vec_true <- t(matrixcalc::elimination.matrix(Rprobit_obj$mod$lRE)) %*% (Rprobit_obj$mod$HO %*% theta_0[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fO)
      } else {
        Omega_vec_est <- (Rprobit_obj$mod$HO %*% Rprobit_obj$theta[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fO)
        Omega_vec_true <- (Rprobit_obj$mod$HO %*% theta_0[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fO)
      }
      Omega_est <- matrix(Omega_vec_est, Rprobit_obj$mod$lRE, Rprobit_obj$mod$lRE)
      Omega_true <- matrix(Omega_vec_true, Rprobit_obj$mod$lRE, Rprobit_obj$mod$lRE)
      out_Omega <- rbind(true = (Omega_true %*% t(Omega_true))[lower.tri(Omega_true, diag = TRUE)], estimate = c(Omega_est %*% t(Omega_est))[lower.tri(Omega_est, diag = TRUE)])
      colnames(out_Omega) <- paste("omega", apply(which(lower.tri(Omega_est, diag = TRUE), arr.ind = T), 1, paste, collapse = ""), sep = "_")
      out <- cbind(out, out_Omega)
    }

    if (Rprobit_obj$mod$lthL > 0) {
      if (Rprobit_obj$mod$ordered == FALSE) {
        nL <- Rprobit_obj$mod$alt
      } else {
        nL <- Rprobit_obj$data_raw$Tp[1]
      }
      if (Rprobit_obj$mod$alt > 1) {
        cur <- Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO
        Sigma_vec_est <- t(elimmat(nL)) %*% (Rprobit_obj$mod$HL %*% Rprobit_obj$theta[(cur + 1):(cur + Rprobit_obj$mod$lthL)] + Rprobit_obj$mod$fL)
        Sigma_vec_true <- t(elimmat(nL)) %*% (Rprobit_obj$mod$HL %*% theta_0[(cur + 1):(cur + Rprobit_obj$mod$lthL)] + Rprobit_obj$mod$fL)
      } else {
        Sigma_vec_est <- (Rprobit_obj$mod$HO %*% Rprobit_obj$theta[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fL)
        Sigma_vec_true <- (Rprobit_obj$mod$HO %*% theta_0[(Rprobit_obj$mod$lthb + 1):(Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO)] + Rprobit_obj$mod$fL)
      }
      Sigma_est <- matrix(Sigma_vec_est, nL, nL)
      Sigma_true <- matrix(Sigma_vec_true, nL, nL)
      out_Sigma <- rbind(true = (Sigma_true %*% t(Sigma_true))[lower.tri(Sigma_true, diag = TRUE)], estimate = (Sigma_est %*% t(Sigma_est))[lower.tri(Sigma_est, diag = TRUE)])
      colnames(out_Sigma) <- paste("sigma", apply(which(lower.tri(Sigma_est, diag = TRUE), arr.ind = T), 1, paste, collapse = ""), sep = "_")
      out <- cbind(out, out_Sigma)
    }
  }

  return(out)
}
