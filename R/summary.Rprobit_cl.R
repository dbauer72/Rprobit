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

  ### build model parameters
  par <- build_par_from_mod(
    theta = object$theta,
    mod = object$mod,
    variances = object$vv
  )

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
    sdall <- sqrt(diag(object$vv))
    sdtauk <- sdall[(length(sdall) - length(par$tauk) + 1):length(sdall)]
    sdtauk[-1] <- sdtauk[-1] * (diff(par$tauk)^2)
    par$tauk_sd <- sdtauk
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
      "prediction_summary" = predict_Rprobit(object),
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
  print_est_sd(
    est_sd = cbind(x$par$Sigma, if(is.null(x$par$Sigma_sd)) x$par$Sigma * NA else x$par$Sigma_sd),
    row_names = x$alt_names,
    col_names = x$alt_names
  )
  if (x$ordered) {
    ### for ordered alternatives: print interval limits
    writeLines("")
    print_est_sd(
      est_sd = matrix(c(x$par$tauk, x$par$tauk_sd), nrow = 1),
      row_names = "tauk",
      col_names = 1:length(x$par$tauk)
    )
  }

  ### Prediction information
  writeLines("")
  print(x$prediction_summary)

  ### Log-likelihood value
  writeLines("")
  writeLines(sprintf("Log-likelihood: %f", x$ll))

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

