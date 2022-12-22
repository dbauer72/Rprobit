#' Initialize fit_Rprobit
#'
#' @description
#' Function that initializes the estimation routine.
#'
#' @param ll_function
#' A \code{function} that computes the log-likelihood value.
#' @param Rprobit_obj
#' An \code{\link{Rprobit_cl}} object.
#' @param init_method
#' A \code{character} specifying the initialization method, one of:
#' * \code{"random"}
#' * \code{"subset"}
#' @param data_tr
#' A transformed data set (utilities are differentiated with respect to the
#' chosen alternative).
#'
#' @return
#' a vector of model parameters
#'
#' @keywords internal

init_Rprobit <- function(ll_function, Rprobit_obj, init_method, data_tr = NULL) {

  ### do not need to calculate hessian for initialization
  Rprobit_obj$control$hess <- 0

  ### count number of initialization trials
  trials <- 0

  ### transform the data
  if (is.null(data_tr)) {
    data_tr <- Rprobit_obj$data$clone()
    data_tr$data <- substract_choice_regressor_from_data(data = Rprobit_obj$data$data)
  }
  ### define parameter dimension
  lth <- Rprobit_obj$mod$lthb + Rprobit_obj$mod$lthO + Rprobit_obj$mod$lthL
  if (Rprobit_obj$mod$ordered == TRUE) {
    lth <- lth + Rprobit_obj$mod$alt - 1 # add parameters for tauk.
  }

  ### get initial theta vector
  if (length(Rprobit_obj$theta) >= lth) {
    theta <- Rprobit_obj$theta[1:lth]
  } else {
    theta <- rep(NA, lth)
  }

  ### check if any initial parameters are NA
  if (!any(!is.na(theta))) theta <- NULL

  ### if an initial theta vector is already given
  if (!is.null(theta)) {
    if (any(is.na(theta))) {
      theta[is.na(theta)] <- stats::rnorm(sum(is.na(theta)))
    }

    if (length(theta) < lth) {
      theta <- stats::rnorm(lth)
      theta[1:length(Rprobit_obj$theta)] <- Rprobit_obj$theta
    }

    if (check_cov_mat(theta = theta, mod = Rprobit_obj$mod)) {
      try_theta <- try(ll_function(theta = theta, data = data_tr, mod = Rprobit_obj$mod, control = Rprobit_obj$control), silent = TRUE)
      if (inherits(try_theta, "try-error")) {
        theta <- NULL
      }
    } else {
      theta <- NULL
    }
  }

  ### if no valid initial theta vector is available
  if (is.null(theta)) {
    while (TRUE) {

      ### draw theta randomly
      if (init_method == "random") {
        theta <- stats::rnorm(lth)
      }

      ### use true theta, if supplied
      if (init_method == "theta") {
        theta <- Rprobit_obj$theta_0
      }

      ### estimate theta on subset of data
      if (init_method == "subset") {
        proportion <- 0.1
        Rprobit_obj$data$data <- Rprobit_obj$data$data[sample.int(n = Rprobit_obj$mod$N, size = round(proportion * length(Rprobit_obj$data)))]
        # Rprobit_obj$mod$N <- length(Rprobit_obj$data$data)
        theta <- fit_Rprobit(Rprobit_obj, init_method = "random")$theta
      }

      ### update counter and test if theta is valid
      trials <- trials + 1
      if (check_cov_mat(theta = theta, mod = Rprobit_obj$mod)) {
        try_theta <- try(ll_function(theta = theta, data = data_tr, mod = Rprobit_obj$mod, control = Rprobit_obj$control), silent = TRUE)
      } else {
        try_theta <- "Error\n"
        attr(try_theta, "class") <- "try-error"
      }

      ### break if initialization was successful or failed too often
      if (trials == 10) {
        stop("Initialization failed.", call. = FALSE)
        break
      }
      if (!inherits(try_theta, "try-error")) {
        break
      }
    }
  }

  return(theta)
}
