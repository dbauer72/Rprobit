#' @rdname Rprobit_cl
#'
#' @param object
#' An \code{\link{Rprobit_cl}} object.
#' @param denormalise
#' A \code{logical}, decides if summary should be given with respect to original
#' variables (\code{denormalise = TRUE}) or normalised variables
#' (\code{denormalise = FALSE}).
#' By default, \code{denormalise = TRUE}.
#' @param ...
#' Not used.
#'
#' @return
#' A \code{summary.Rprobit_cl} object, which is a \code{list} that contains
#' the following information on the estimated model:
#' * \code{N}, the number of decision makers
#' * \code{alt}, the number of alternatives
#' * \code{Tp}, the vector of number of choice occasions
#' * \code{omega_present}, a \code{logical}, indicating whether an Omega matrix exists
#' * \code{vars}, a \code{character} vector of variable names
#' * \code{alt_names}, a \code{character} vector of alternative names
#' * \code{ordered}, a \code{logical}, indicating whether the alternatives are ordered
#' * \code{par}, a \code{list} of parameter estimates
#' * \code{prediction_summary}, a \code{matrix} of prediction details
#' * \code{ll}, a \code{numeric}, the log-likelihood value
#'
#' @exportS3Method

summary.Rprobit_cl <- function(object, denormalise = TRUE, ...) {

  LC_model = FALSE
  StSp_model = FALSE
  ### build model parameters
  if (!is.null(object$mod$num_class)){
    LC_model = TRUE
    num_class <- object$mod$num_class
    par_all <- list(num_class)
    mod <- object$mod
    for (j in 1:num_class){
        param_one <- mod$lthb + mod$lthO + mod$lthL 
        ind_j <- (j-1)*param_one+c(1:param_one)
        par_all[[j]] <- build_par_from_mod(object$theta[ind_j], mod, variances = object$vv[ind_j,ind_j,drop=FALSE])
    }
    pi_eta <- c(0,object$theta[num_class*param_one + c(1:(num_class-1))])
    pi <- exp(pi_eta)/sum(exp(pi_eta))
    ind_gam <- num_class*param_one + c(1:(num_class-1)) 
    pi_V <- object$vv[ind_gam,ind_gam]
    
    par <- list( par_all = par_all,pi = pi, pi_V =  pi_V)
  } else {
    par <- build_par_from_mod(
      theta = object$theta,
      mod = object$mod,
      variances = object$vv
    )
  }
  ### check if Omega matrix is available
  omega_present <- FALSE
  if (!any(is.na(object$mod$HO))) {
    if (is.matrix(object$mod$HO)) {
      if (nrow(object$mod$HO) > 0) {
        omega_present <- TRUE
      }
    }
  }

  ### calculate back from normalised data
  denormalise_check <- is.list(object$control$normalisation) &
    (length(grep("ASC_", object$vars)) > 0) & denormalise

  if (denormalise_check) {

    ### get the ASC of the reference alternative implied by the normalisation
    ref_alt <- which(!(paste("ASC_", object$alt_names, sep = "") %in% object$vars))
    non_ref_alt <- which((paste("ASC_", object$alt_names, sep = "") %in% object$vars))
    A_matrix <- diag(length(object$vars))

    for (i_var in object$vars) {
      norm_list <- object$control$normalisation[[i_var]]
      if (is.null(norm_list)) {
        norm_list <- object$control$normalisation[[paste0(i_var, "_", object$alt_names[ref_alt])]]
      }
      if (is.null(norm_list)) {
        norm_names <- names(object$control$normalisation)
        norm_var_which <- which(paste0(rep(norm_names, length(object$alt_names)), "_", rep(object$alt_names, each = length(norm_names)), sep = "") == i_var)
        if (length(norm_var_which) > 0) {
          norm_name <- rep(norm_names, length(object$alt_names))[norm_var_which[1]]
          norm_list <- object$control$normalisation[[norm_name]]
        }
      }
      if (!is.null(norm_list)) {
        a_norm <- norm_list$a_norm
        b_norm <- norm_list$b_norm
        A_matrix[which(object$vars == i_var), which(object$vars == i_var)] <- 1 / b_norm
        for (i_alt in 1:length(object$alt_names)) {
          alt_connection <- grep(paste0("_", object$alt_names)[i_alt], i_var)
          if (length(alt_connection) > 0) {
            if (i_alt == ref_alt) {
              A_matrix[object$vars %in% paste("ASC_", object$alt_names[non_ref_alt], sep = ""), which(object$vars == i_var)] <- A_matrix[object$vars %in% paste("ASC_", object$alt_names[non_ref_alt], sep = ""), which(object$vars == i_var)] + a_norm / b_norm
            } else {
              A_matrix[object$vars == paste("ASC_", object$alt_names[i_alt], sep = ""), which(object$vars == i_var)] <- A_matrix[object$vars == paste("ASC_", object$alt_names[i_alt], sep = ""), which(object$vars == i_var)] - a_norm / b_norm
            }
          }
        }
      }
    }

    ### Calculate variance of the estimates, considering covariance between the regressors
    if (object$control$hess == 1) {
      b_cov_mat <- object$mod$Hb %*% object$vv[1:object$mod$lthb, 1:object$mod$lthb] %*% t(object$mod$Hb)
      par$b_sd <- sqrt(diag(A_matrix %*% b_cov_mat %*% t(A_matrix)))
    }

    ### TODO: Omega_sd_print does not yet incorporate covariances between the parameters theta_Omega
    par$Omega <- A_matrix %*% par$Omega %*% t(A_matrix)
    if (!is.null(par$Omega_sd)) {
      par$Omega_sd <- A_matrix %*% par$Omega_sd %*% t(A_matrix)
    }
  }

  if (object$mod$ordered) {
    # variance for treshholds
    ind_thresh <- (object$mod$lthb +object$mod$lthO + object$mod$lthL + c(1:object$mod$alt-1))
    l_it <- length(ind_thresh)
    Vthresdh <- object$vv[ind_thresh,ind_thresh]
    
    dtauk <- diff(par$tauk)
    L <- matrix(0,l_it,l_it)
    L[,1] <- 1
    for (j in 2:l_it){
      L[j:l_it,j] <- dtauk[j-1]
    }
    par$tauk_sd <- sqrt(diag( L %*% Vthresdh %*% t(L)))
    
    # for state space models include matrices. 
    if (inherits(object$mod,"mod_StSp_cl")){
      n = object$mod$dim_state
      s = dim(par$Sigma)[1]
      ind_StSp <- (object$mod$lthb +object$mod$lthO + object$mod$lthL + object$mod$alt-1 + c(1:(2*s*n)))

      param <- object$theta[ind_StSp]
      sd_param <- sqrt(diag(object$vv[ind_StSp,ind_StSp]))
      stsp <- param_to_system_R(param, s=s, n=n, grad_bool=0, stationary = FALSE)
      sd_stsp <- param_to_system_R(sd_param, s=s, n=n, grad_bool=1, stationary = FALSE)
      
      # add to output 
      par$StSp <- stsp
      par$sd_StSp <- sd_stsp
      
      StSp_model = TRUE
    }
    # for AR models include matrices
    # for state space models include matrices. 
    if (inherits(object$mod,"mod_AR_cl")){
      n = object$mod$lag_length
      ind_StSp <- (object$mod$lthb +object$mod$lthO + object$mod$lthL + object$mod$alt-1 + c(1:n))
      
      param <- object$theta[ind_StSp]
      sd_param <- sqrt(diag(object$vv[ind_StSp,ind_StSp]))
      stsp <- param_to_system_AR_R(param, lag=n, grad_bool=0, stationary = FALSE)
      sd_stsp <- param_to_system_AR_R(sd_param, lag=n, grad_bool=1, stationary = FALSE)
      
      # add to output 
      par$StSp <- stsp
      par$sd_StSp <- sd_stsp
      
      StSp_model = TRUE
    }
    
  }

  ### build output
  structure(
    list(
      "N" = object$data_raw$N,
      "alt" = object$mod$alt,
      "Tp" = object$data_raw$Tp,
      "omega_present" = omega_present,
      "vars" = object$vars,
      "alt_names" = object$alt_names,
      "ordered" = object$mod$ordered,
      "par" = par,
      "latent_class_model" = LC_model,
      "state_space_model" = StSp_model,
      "prediction_summary" = predict_Rprobit(object),
      "predictions" = predict_Rprobit(object,all_pred =TRUE),
      "ll" = object$ll
    ),
    class = c("summary.Rprobit_cl", "list")
  )
}

#' @noRd
#' @exportS3Method

print.summary.Rprobit_cl <- function(x, ...) {

  writeLines("Summary of probit model:")

  ### Data information
  writeLines("")
  writeLines(paste("N =", x$N, "decision makers"))
  writeLines(paste("alt =", x$alt, "alternatives"))
  if (length(unique(x$Tp)) == 1) {
    writeLines(paste("Tp =", unique(x$Tp), "choice occasions"))
  }
  if (length(unique(x$Tp)) > 1) {
    writeLines(paste("Tp =", min(x$Tp), "to", max(x$Tp), "choice occasions"))
  }
  writeLines("")

  ### Model information
  if (x$latent_class_model){
    num_class <- length(x$par[[1]])
    writeLines(paste("Latent class model with ",length(x$par[[1]])," classes."))
    
    writeLines("Class membership percentages:")
    writeLines("pi =")
    pi <- matrix(x$par[[2]],ncol=1)
    # calculate variance of estimated percentages 
    dpidgamma <- diag(as.vector(pi)) - pi %*% t(pi)
    dpidgamma <- dpidgamma[,-1, drop = FALSE] 
    V_pi <- dpidgamma %*% x$par$pi_V %*% t(dpidgamma)
    
    est_sd = cbind(t(pi), t(sqrt(diag(V_pi))))
    
    print_est_sd(
      est_sd,
      row_names = "Perc",
      col_names = paste("Cl. ",1:num_class)
    )
    writeLines("")
  }
  writeLines("U = X * beta + eps")
  if (x$omega_present) {
    writeLines("beta ~ MVN(b,Omega)")
  } else {
    writeLines("beta = b")
  }
  writeLines("eps ~ MVN(0,Sigma)")
  writeLines("")
  writeLines("b =")

  ### Parameter information
  if (x$latent_class_model == FALSE){ 
    # only one class, standard MNP model 
    print_est_sd(
      est_sd = cbind(x$par$b, if(is.null(x$par$b_sd)) x$par$b * NA else x$par$b_sd),
      row_names = x$vars
    )
    writeLines("")
    if (x$omega_present) {
      writeLines("Omega =")
      print_est_sd(
        est_sd = cbind(x$par$Omega, if(is.null(x$par$Omega_sd)) x$par$Omega * NA else x$par$Omega_sd),
        row_names = x$vars,
        col_names = x$vars
      )
      writeLines("")
    }
    writeLines("Sigma =")
    if (x$ordered == FALSE) { 
      print_est_sd(
        est_sd = cbind(x$par$Sigma, if(is.null(x$par$Sigma_sd)) x$par$Sigma * NA else x$par$Sigma_sd),
        row_names = x$alt_names,
        col_names = x$alt_names
      )
    }
  } else {
    par_all <- x$par[[1]]
    b_mat <- par_all[[1]]$b
    for (j in 2:num_class){
      b_mat <- cbind(b_mat,par_all[[j]]$b)
    }
    
    b_mat <- cbind(b_mat,if(is.null(par_all[[1]]$b_sd)) par_all[[1]]$b * NA else par_all[[1]]$b_sd )
    for (j in 2:num_class){
      b_mat <- cbind(b_mat, if(is.null(par_all[[j]]$b_sd)) par_all[[j]]$b * NA else par_all[[j]]$b_sd )
    }
    
    print_est_sd(
      est_sd = cbind(b_mat),
      row_names = x$vars,
      col_names = paste("Cl. ",1:num_class)
    )
    
    writeLines("")
    if (x$omega_present) {
      writeLines("Omega =")
      print_est_sd(
        est_sd = cbind(par_all[[1]]$Omega, if(is.null(par_all[[1]]$Omega_sd)) par_all[[1]]$Omega * NA else par_all[[1]]$Omega_sd),
        row_names = x$vars,
        col_names = x$vars
      )
      writeLines("")
    }
    writeLines("Sigma =")
    print_est_sd(
      est_sd = cbind(par_all[[1]]$Sigma, if(is.null(par_all[[1]]$Sigma_sd)) par_all[[1]]$Sigma * NA else par_all[[1]]$Sigma_sd),
      row_names = x$alt_names,
      col_names = x$alt_names
    )
  }
  
  if (x$ordered) {
    Tp <- dim(x$par$Sigma)[1]
    print_est_sd(
      est_sd = cbind(x$par$Sigma, if(is.null(x$par$Sigma_sd)) x$par$Sigma * NA else x$par$Sigma_sd),
      row_names = 1:Tp,
      col_names = 1:Tp
    )
    ### for ordered alternatives: print interval limits
    writeLines("")
    print_est_sd(
      est_sd = matrix(c(x$par$tauk, x$par$tauk_sd), nrow = 1),
      row_names = "tauk",
      col_names = 1:length(x$par$tauk)
    )
    
    if (x$state_space_model){
      n = dim(x$par$StSp$A)[1]
      s = dim(x$par$StSp$C)[1]
      
      writeLines("State Space Model for Error Autocorrelation")
      writeLines("A =")      
      print_est_sd(
        est_sd = cbind(x$par$StSp$A, if(is.null(x$par$sd_StSp$A)) x$par$StSp$A * NA else x$par$sd_StSp$A),
        row_names = 1:n,
        col_names = 1:n
      )
      writeLines("K =")      
      print_est_sd(
        est_sd = cbind(x$par$StSp$K, if(is.null(x$par$sd_StSp$K)) x$par$StSp$K * NA else x$par$sd_StSp$K),
        row_names = 1:n,
        col_names = 1:s
      )
      writeLines("C =")      
      print_est_sd(
        est_sd = cbind(x$par$StSp$C, if(is.null(x$par$sd_StSp$C)) x$par$StSp$C * NA else x$par$sd_StSp$C),
        row_names = 1:s,
        col_names = 1:n
      )
    }
      
    
  }
  ### Prediction information
  writeLines("")
  print(x$prediction_summary)

  ### Log-likelihood value
  writeLines("")
  writeLines(sprintf("Log-CML value: %f", x$ll))

}

#' Print Estimated and Standard Deviations
#'
#' @description
#' This helper function prints a formatted version of a matrix with estimates
#' and standard deviations.
#'
#' @param est_sd
#' A \code{matrix} of dimension \code{m} times \code{2n}, where the first
#' \code{m} columns are estimates and the
#' last \code{m} columns are corresponding standard deviations.
#' @param row_names
#' A \code{character} vector of length \code{m} with row names.
#' Can be \code{NULL} (default) for no row names.
#' @param col_names
#' A \code{character} vector of length \code{n} with column names.
#' Can be \code{NULL} (default) for no column names.
#'
#' @return
#' No return value, prints a \code{matrix}.
#'
#' @examples
#' \dontrun{
#' est_sd <- matrix(1:12, nrow = 3, ncol = 4)
#' print_est_sd(est_sd, row_names = LETTERS[1:3], col_names = LETTERS[4:5])
#' }
#'
#' @keywords internal utils

print_est_sd <- function(est_sd, row_names = NULL, col_names = NULL) {

  ### compute dimensions
  nr <- dim(est_sd)[1]
  nc <- dim(est_sd)[2] / 2

  ### format the rownames
  if (is.null(row_names)) {
    row_names_print <- rep("", nr)
  } else {
    row_names_print <- row_names
    max_length <- max(nchar(row_names_print))
    for (i in 1:length(row_names_print)) {
      while (nchar(row_names_print[i]) < max_length) {
        row_names_print[i] <- paste(row_names_print[i], " ", sep = "")
      }
    }
  }

  ### format the entries to have equal length
  est_sd_print <- est_sd
  for (r in 1:nr) {
    for (c in 1:(2 * nc)) {
      est_sd_print[r, c] <- sprintf("%f", est_sd[r, c])
    }
  }
  for (c in 1:(2 * nc)) {
    max_length <- max(nchar(est_sd_print[, c]))
    for (r in 1:nr) {
      while (nchar(est_sd_print[r, c]) < max_length) {
        est_sd_print[r, c] <- paste(" ", est_sd_print[r, c], sep = "")
      }
    }
  }

  ### add column names
  if (!is.null(col_names)) {
    col_names_print <- as.character(col_names)
    str <- paste(rep(" ", nchar(row_names_print[1]) + 1), sep = "", collapse = "")
    for (c in 1:nc) {
      while (nchar(col_names_print[c]) < (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 4)) {
        col_names_print[c] <- paste(col_names_print[c], " ", sep = "")
      }
      if (nchar(col_names_print[c]) > (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 5)) {
        col_names_print[c] <- paste(substring(col_names_print[c], 1, (nchar(est_sd_print[1, c]) + nchar(est_sd_print[1, c + nc]) + 2)), ". ", sep = "")
      }

      str <- c(str, col_names_print[c])
    }
    cat(c(str, "\n"))
  }

  ### print the matrix row by row
  for (r in 1:nr) {
    str <- row_names_print[r]
    for (c in 1:nc) {
      str <- c(str, sprintf(" %s (%s)", est_sd_print[r, c], est_sd_print[r, c + nc]))
    }
    cat(c(str, "\n"))
  }
}

