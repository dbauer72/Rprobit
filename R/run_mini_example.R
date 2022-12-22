#' Perform mini simulation examples
#'
#' @description This function performs mini simulation examples.
#'
#' @param form
#' A \code{form} object.
#' @param mod
#' A \code{\link{mod_cl}} object.
#' @param re
#' A \code{character} vector of variable names which are considered to have
#' random effects.
#' @param approx_method
#' A \code{character}. Can be (one or more of)
#' * \code{"SJ"} (default),
#' * \code{"TVBS"},
#' * \code{"TVBSv2"},
#' * \code{"ME"}.
#' @param show_theta boolean, if true, compare theta vector, if false (default) compare entries of b, Omega and Sigma
#' @param seed set a seed for the initial values and the simulated data, no seed per default
#' @param at_true if \code{TRUE}, starts optimization from true parameter
#' @param runs number of numerical optimizations at different starting values, default \code{1}
#' @param probit if \code{TRUE} (default) estimates a probit model, else a CML-model
#' @param nCores number of cores used for estimation
#' @param print.level Passed on to \code{nlm}, default \code{0}
#' @param normalize boolean indicating whether probabilities are normalize to sum to 1, passed on to \code{nlm}, default \code{FALSE}
#' @param control_simul list containing at least field Tp controlling the simulation.
#' @param verbose boolean indicating whether to print on-screen information
#' @param compare_time boolean indicating whether to save estimation time
#' @param theta parameter vector
#' @param show_ll boolean indicating whether to save the log-likelihood
#' @param cml_pair_type integer specifying the CML to use (0: full pairwise, 1: adjacent lopped, 2: adjacent)
#' By default \code{cml_pair_type = 0} is used
#' @return no return value; on-screen information
#'
#' @importFrom stats rnorm runif
#' 
#' @keywords internal

run_mini_example <- function(
    form, mod, re, approx_method = "SJ", show_theta = FALSE, seed = NULL,
    at_true = FALSE, runs = 1, probit = TRUE, nCores = 1, print.level = 0,
    normalize = FALSE, control_simul = NULL, verbose = TRUE, compare_time = FALSE, theta = NULL,
    show_ll = FALSE, cml_pair_type = 0
 ) {

  ### check input 'approx_method'
  if (!(is.character(approx_method) &&
        all(approx_method %in% c("SJ", "TVBS", "TVBSv2", "ME")))) {
    stop("Input 'approx_method' misspecified.", call. = FALSE)
  }

  ### draw true model parameters
  if (!is.null(seed)) {
    set.seed(seed)
  }

  ### draw true parameter vector
  if (is.null(theta)){
    theta_0 <- if (mod$ordered) {
      c(stats::rnorm(mod$lthb + mod$lthO),
        stats::runif(mod$lthL, 1, 2),
        -2,
        stats::runif(mod$alt - 2, -1, 1)
        )
    } else {
      stats::rnorm(mod$lthb + mod$lthO + mod$lthL)
    }
  } else {
    theta_0 <- theta
  }

  ### create 'Rprobit_obj' object
  Rprobit_obj <- setup_Rprobit(
    form = form,
    mod = mod,
    re = re,
    seed = seed,
    theta_0 = theta_0,
    control = control_simul
  )
  Rprobit_obj$control$hess <- 1
  Rprobit_obj$control$probit <- probit
  Rprobit_obj$control$nCores <- nCores
  Rprobit_obj$control$normalize <- normalize

  # add CML type
  Rprobit_obj$control$control_weights$cml_pair_type = cml_pair_type 
  
  ### save estimates
  estimates <- list()

  ### loop over runs
  for (run in 1:runs) {

    if (verbose) {
      writeLines(paste0("Optimization run ", run, ":"))
    }

    ### specify initial values
    if (at_true) {
      Rprobit_obj$theta <- theta_0
    } else {
      Rprobit_obj$theta <- init_Rprobit(ll_function = ll_macml, Rprobit_obj = Rprobit_obj, init_method = "random")
    }

    ### optimize wrt "approx_method"
    comp <- matrix(theta_0, nrow = 1)
    if(compare_time){
      time_est <- numeric(length(approx_method)+1)
      Rprobit_obj$info$estimation_time <- NULL
      Rprobit_obj$info$hess_time <- NULL
    }
    if(show_ll){
      ll_est <- numeric(length(approx_method)+1)
      Rprobit_obj$ll <- NA
    }
    rownames(comp) <- "true"
    for (am in approx_method) {
      Rprobit_obj$control$approx_method <- am
      if (at_true) {
        Rprobit_obj$theta <- theta_0
        Rprobit_obj <- fit_Rprobit(Rprobit_obj = Rprobit_obj, init_method = "theta")
      } else {
        Rprobit_obj <- fit_Rprobit(Rprobit_obj = Rprobit_obj)
      }
      if(compare_time){
        time_est[which(approx_method == am)+1] <- Rprobit_obj$info$estimation_time + Rprobit_obj$info$hess_time
      }
      if(show_ll){
        ll_est[which(approx_method == am)+1] <- Rprobit_obj$ll
      }
      comp <- rbind(comp, Rprobit_obj$theta)
      rownames(comp)[nrow(comp)] <- am
      if (verbose) {
        cat("Approximation method:", am, "\n")
        print(compare_sim_fit(theta_0, Rprobit_obj, show_theta = show_theta))
        cat("\n")
      }
    }

    ### add mean squared error
    comp <- cbind(comp, "mse" = apply(comp, 1, function(x) mean((x - comp[1, ])^2)))

    ### add time
    if(compare_time){
      comp <- cbind(comp, "time" = time_est)
      rm(time_est)
    }
    ### add ll
    if(show_ll){
      comp <- cbind(comp, "ll" = ll_est)
      rm(ll_est)
    }
    
    estimates[[run]] <- comp

    ### print summary of fitted model
    if (verbose) {
      summary(Rprobit_obj)
    }
  }

  return(estimates)
}
