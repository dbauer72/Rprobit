#' Calculate unbalanced-panel weights
#' @description Function that calculates unbalanced-panel weights from Rprobit_obj object
#' @param Rprobit_obj \code{Rprobit_obj}-object
#' @return dataframe with weights depending on number of observations
#' @keywords internal

build_unbalanced_panel_weights <- function(Rprobit_obj) {
  Tp <- Rprobit_obj$data_raw$Tp

  ### check if data is actually unbalanced panel data
  if (length(unique(Tp)) == 1) {
    print("Data resembles a balanced panel")

    ### get number of pairs per individual
    if (!is.null(Rprobit_obj$control$pairs_list)) {
      C_s <- nrow(Rprobit_obj$control$pairs_list[[1]])
    } else if (!is.null(Rprobit_obj$control$control_weights$cml_pair_type)) {
      if (Rprobit_obj$control$control_weights$cml_pair_type == 0) {
        C_s <- Tp[1] * (Tp[1] - 1) / 2
      } else if (Rprobit_obj$control$control_weights$cml_pair_type == 1) {
        C_s <- Tp[1]
      } else if (Rprobit_obj$control$control_weights$cml_pair_type == 2) {
        C_s <- Tp[1] - 1
      } 
    } else if (!is.null(Rprobit_obj$control$cml_pair_type)) {
      if (Rprobit_obj$control$cml_pair_type == 0) {
        C_s <- Tp[1] * (Tp[1] - 1) / 2
      } else if (Rprobit_obj$control$cml_pair_type == 1) {
        C_s <- Tp[1]
      } else if (Rprobit_obj$control$cml_pair_type == 2) {
        C_s <- Tp[1] - 1
      } 
    } else {
      C_s <- Tp[1]
    }

    ### return weight equal to one
    out <- data.frame(n_obs = Tp[1], weights = 1, n_ind = Rprobit_obj$mod$N, n_pairs = C_s)
  } else {
    ### if data has unbalanced panel structure

    ### extract possible number of observations
    Tp_s <- sort(unique(Tp))

    ### initialize vectors
    s <- numeric(length(Tp_s))
    C_s <- s
    N_s <- s
    v_s <- s

    ### go through possible numbers of observations
    for (i_tp in 1:length(Tp_s)) {

      ### get number of observations
      s[i_tp] <- Tp_s[i_tp]

      ### get individuals involved
      n_s <- which(Rprobit_obj$data_raw$Tp == s[i_tp])

      ### get number of individuals
      N_s[i_tp] <- length(n_s)

      ### get number of pairs per individual
      if (!is.null(Rprobit_obj$control$pairs_list)) {
        C_s[i_tp] <- nrow(Rprobit_obj$control$pairs_list[[n_s[1]]])
      } else if (!is.null(Rprobit_obj$control$control_weights$cml_pair_type)) {
        if (Rprobit_obj$control$control_weights$cml_pair_type == 0) {
          C_s[i_tp] <- s[i_tp] * (s[i_tp] - 1) / 2
        } else if (Rprobit_obj$control$control_weights$cml_pair_type == 1) {
          C_s[i_tp] <- s[i_tp]
        } else if (Rprobit_obj$control$control_weights$cml_pair_type == 2) {
          C_s[i_tp] <- s[i_tp] - 1
        }
      } else {
        C_s[i_tp] <- s[i_tp]
      }
    }


    ### determine which weighting method should be used
    method <- "BB"
    if (!is.null(Rprobit_obj$control$control_weights$unbalanced_panel_weights$method)) {
      method <- Rprobit_obj$control$control_weights$unbalanced_panel_weights$method
    }

    ### apply weighting method
    if (method %in% c("JoeLee", "JL")) {
      rho <- 0.5
      if (!is.null(Rprobit_obj$control$control_weights$unbalanced_panel_weights$rho)) {
        rho <- Rprobit_obj$control$control_weights$unbalanced_panel_weights$rho
      }

      ### calculate JoeLee weights
      w_s <- 1 / ((s - 1) * (1 + (s - 1) * rho))
      sum_weight <- sum(w_s * C_s * N_s)
      target_weight <- sum(s * N_s)
      w_s_scaled <- w_s / sum_weight * target_weight

      ### prepare output
      out <- data.frame(n_obs = s, weights = w_s_scaled, n_ind = N_s, n_pairs = C_s)
    } else {
      ### use BB weights

      ### extract variance-metric function
      var_metric <- function(Vs) {
        sum(diag(Vs))
      }
      if (!is.null(Rprobit_obj$control$control_weights$unbalanced_panel_weights$var_metric)) {
        var_metric <- Rprobit_obj$control$control_weights$unbalanced_panel_weights$var_metric
      }

      ### check if Hessian matrix exists in Rprobit_obj
      if (any(is.na(Rprobit_obj$H))) {
        ### calculate hessian (or approximate hessian)
        # print("Estimate hessian matrix H")

        ### check if hessian or approximate hessian should be used
        hess <- FALSE
        if (!is.null(Rprobit_obj$control$control_weights$unbalanced_panel_weights$approx_hess)) {
          hess <- 1 - as.numeric(Rprobit_obj$control$control_weights$unbalanced_panel_weights$approx_hess)
        }
        Rprobit_obj$control$hess <- hess
        Rprobit_obj$control$el <- TRUE
        save_control_weights <- Rprobit_obj$control$control_weights$unbalanced_panel_weights
        Rprobit_obj$control$control_weights$unbalanced_panel_weights <- NULL
        Rprobit_obj <- fit_Rprobit(Rprobit_obj = Rprobit_obj)
        Rprobit_obj$control$control_weights$unbalanced_panel_weights <- save_control_weights
      }
      if (is.null(Rprobit_obj$gradEL)) {
        ### calculate gradient gradient contribution of each individual
        print("Calculate gradient gradient contribution of each individual")
        Rprobit_obj$control$el <- TRUE
        save_control_weights <- Rprobit_obj$control$control_weights$unbalanced_panel_weights
        Rprobit_obj$control$control_weights$unbalanced_panel_weights <- NULL
        Rprobit_obj$control$control_nlm$force_hess <- TRUE
        Rprobit_obj <- fit_Rprobit(Rprobit_obj = Rprobit_obj)
        Rprobit_obj$control$control_weights$unbalanced_panel_weights <- save_control_weights
      }

      ### calculate matrix H_0 and its inverse
      H_0 <- Rprobit_obj$H / sum(N_s * C_s)
      H_0_inv <- adjust_variance_matrix(H_0)$mat_inv
      # H_0 <- adjust_variance_matrix(H_0)$mat
      # H_0_inv <- solve(H_0)

      ### go again through possible numbers of observations
      for (i_tp in 1:length(Tp_s)) {
        ### get individuals involved
        n_s <- which(Rprobit_obj$data_raw$Tp == s[i_tp])

        ### calculate V_s
        V_s <- Rprobit_obj$H * 0
        for (i_ns in n_s) {
          V_s <- V_s + Rprobit_obj$gradEL[[i_ns]] %*% t(Rprobit_obj$gradEL[[i_ns]])
        }

        ### extract v_s
        v_s[i_tp] <- var_metric(H_0_inv %*% V_s %*% H_0_inv)
      }

      ### calculate optimal weights
      w_s <- C_s * N_s / v_s
      sum_weight <- sum(w_s * C_s * N_s)
      target_weight <- sum(s * N_s)
      w_s_scaled <- w_s / sum_weight * target_weight

      ### prepare output
      out <- data.frame(n_obs = s, weights = w_s_scaled, n_ind = N_s, n_pairs = C_s)
      attr(x = out, which = "theta") <- Rprobit_obj$theta
    }
  }

  return(out)
}
