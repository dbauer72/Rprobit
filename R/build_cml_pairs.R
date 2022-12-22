#' Build CML Weighting Matrices
#'
#' @description
#' This function builds a list of matrices to define CML pairs and weights.
#'
#' @param Rprobit_obj
#' An \code{\link{Rprobit_cl}} object.
#'
#' @param cml_pair_type
#' A named \code{list} of parameters defining CML pair weights with
#' - range parameter \code{"d"},
#' - \code{list} of positions of observations \code{"positions"},
#' - \code{"weighting_function"} or weighting strategy \code{"name"},
#' - individual_weights as a \code{vector} of individual weights.
#' By default, \code{cml_pair_type = NULL}.
#'
#' @keywords internal
#'
#' @return
#' A \code{list} of matrices defining CML pairs and weights.

build_cml_pairs <- function(Rprobit_obj, cml_pair_type = NULL) {

  ### get number of observations for each individual
  Tp <- Rprobit_obj$data_raw$Tp

  ### initialize output
  out <- list()

  if (is.null(cml_pair_type)) {
    if (!is.null(Rprobit_obj$control$control_weights$cml_pair_type)) {
      cml_pair_type <- Rprobit_obj$control$control_weights$cml_pair_type
    } else {
      cml_pair_type <- 0
    }
  }

  if (!is.list(cml_pair_type)) {
    if (is.numeric(cml_pair_type) | is.character(cml_pair_type)) {
      cml_pair_type <- list("name" = cml_pair_type)
    }
  }

  if (is.list(cml_pair_type)) {

    ### get distance parameter
    d <- cml_pair_type$d
    ### extract potential list of positions
    positions <- cml_pair_type$positions
    ### extract weighting function
    weighting_function <- cml_pair_type$weighting_function
    ### extract individual weights
    individual_weights <- cml_pair_type$individual_weights
    ### extract unbalanced panel weights
    unbalanced_panel_weights <- cml_pair_type$unbalanced_panel_weights

    ### if no weighting function is defined, check if a name for a cml type was given
    if (is.null(weighting_function)) {
      if (!is.null(cml_pair_type$name)) {
        ### define different decay/growth functions for given names

        if (cml_pair_type$name == "step_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(as.numeric(abs(positions[p_1] - positions[p_2]) <= d))
          }
        } else if (cml_pair_type$name == "step_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - as.numeric(abs(positions[p_1] - positions[p_2]) <= d))
          }
        } else if (cml_pair_type$name == "triangular_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(as.numeric(abs(positions[p_1] - positions[p_2]) <= d) * (1 - abs(positions[p_1] - positions[p_2]) / d))
          }
        } else if (cml_pair_type$name == "triangular_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - as.numeric(abs(positions[p_1] - positions[p_2]) < d) * (1 - abs(positions[p_1] - positions[p_2]) / d))
          }
        } else if (cml_pair_type$name == "triangular_growth_gap") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return((1 - as.numeric(abs(positions[p_1] - positions[p_2]) <= d)) * (1 - as.numeric(abs(positions[p_1] - positions[p_2] + d) <= d) * (1 - abs(positions[p_1] - positions[p_2] + d) / d)))
          }
        } else if (cml_pair_type$name == "Epanechnikov_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(as.numeric(abs(positions[p_1] - positions[p_2]) <= d) * (1 - abs(positions[p_1] - positions[p_2]) / d)^2)
          }
        } else if (cml_pair_type$name == "Epanechnikov_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - as.numeric(abs(positions[p_1] - positions[p_2]) <= d) * (1 - abs(positions[p_1] - positions[p_2]) / d)^2)
          }
        } else if (cml_pair_type$name == "Epanechnikov_growth_gap") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return((1 - as.numeric(abs(positions[p_1] - positions[p_2]) <= d)) * (1 - as.numeric(abs(positions[p_1] - positions[p_2] + d) <= d) * (1 - abs(positions[p_1] - positions[p_2] + d) / d)^2))
          }
        } else if (cml_pair_type$name == "exponential_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(exp(-(abs(positions[p_1] - positions[p_2]) / d)^1 * log(2)))
          }
        } else if (cml_pair_type$name == "exponential_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - exp(-(abs(positions[p_1] - positions[p_2]) / d)^1 * log(2)))
          }
        } else if (cml_pair_type$name == "Weibull_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(exp(-(abs(positions[p_1] - positions[p_2]) / d)^2 * log(2)))
          }
        } else if (cml_pair_type$name == "Weibull_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - exp(-(abs(positions[p_1] - positions[p_2]) / d)^2 * log(2)))
          }
        } else if (cml_pair_type$name == "Hill_decay") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 / (1 + (abs(positions[p_1] - positions[p_2]) / d)^2))
          }
        } else if (cml_pair_type$name == "Hill_growth") {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(1 - 1 / (1 + (abs(positions[p_1] - positions[p_2]) / d)^2))
          }
        } else if (cml_pair_type$name %in% c(1, "adjacent_looped")) {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(as.numeric((p_2 == p_1 + 1) | (p_1 == 1 & p_2 == length(positions))))
          }
        } else if (cml_pair_type$name %in% c(2, "adjacent_chained", "adjacent")) {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(as.numeric(p_2 == p_1 + 1))
          }
        } else {
          weighting_function <- function(p_1, p_2, d, positions) {
            return(rep(1, length(p_1)))
          }
        }
      } else {
        weighting_function <- function(p_1, p_2, d, positions) {
          return(rep(1, length(p_1)))
        }
      }
    }



    ### build pairs matrices
    for (i_n in 1:Rprobit_obj$data_raw$N) {
      Tp_n <- Rprobit_obj$data_raw$Tp[i_n]
      pairs_matrix_n <- matrix(0, Tp_n * (Tp_n - 1) / 2, 3)
      if (Tp_n>1){
  
        if (is.list(positions)) {
          if (is.numeric(positions[[i_n]])) {
            positions_n <- positions[[i_n]]
          }
        } else {
          positions_n <- 1:Tp_n
        }

        i_pair <- 1
        for (i_1 in 1:(Tp_n - 1)) {
          for (i_2 in (i_1 + 1):Tp_n) {
            pairs_matrix_n[i_pair, 1:2] <- c(i_1, i_2) - 1
            pairs_matrix_n[i_pair, 3] <- weighting_function(i_1, i_2, d, positions_n)

            i_pair <- i_pair + 1
          }
        }
        pairs_matrix_n <- pairs_matrix_n[pairs_matrix_n[, 3] > 0, , drop = FALSE]

        ### if individual weights are defined, multiply weights by individual weights
        if (!is.null(individual_weights)) {
          pairs_matrix_n[, 3] <- pairs_matrix_n[, 3] * individual_weights[i_n]
        }

        ### if unbalanced panel weights are defined, multiply weights by individual weights
        if (!is.null(unbalanced_panel_weights)) {
          pairs_matrix_n[, 3] <- pairs_matrix_n[, 3] * unbalanced_panel_weights[i_n]
        }
      }
      

      out[[i_n]] <- pairs_matrix_n
    } ### end for i_n
  } ### end if is.list(cml_pair_type)


  ### if no list was successfully built, use void type
  if (length(out) == 0) {
    ### when CML pairs are a standard type
    out["NP_LIST"] <- 0
  }


  return(out)
}
