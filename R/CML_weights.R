#' Calculate and Apply CML Weights
#'
#' @description
#' This function calculates and applies CML weights for probit modeling.
#'
#' @param Rprobit_obj
#' A \code{\link{Rprobit_cl}} object.
#' @param control_weights
#' A named \code{list} of control parameters to specify the CML weightings.
#'
#' @return
#' A \code{\link{Rprobit_cl}} object.
#'
#' @export

CML_weights <- function(Rprobit_obj, control_weights = NULL) {
  ### if no specific control weights are supplied, use those from Rprobit_obj
  if (is.null(control_weights)) {
    control_weights <- Rprobit_obj$control$control_weights
  }

  ### check if any weight controls are specified or if pairs_list is already defined
  individual_weights <- NULL
  if (is.null(control_weights) | !is.null(Rprobit_obj$control$pairs_list)) {
    ### TODO: should control_weights setting be set to NULL if pairs_list is used instead or not?
    Rprobit_obj$control$control_weights <- NULL
    out <- Rprobit_obj
  } else {

    unbalanced_panel_weights <- Rprobit_obj$data_raw$Tp/max(Rprobit_obj$data_raw$Tp)
    control_weights$unbalanced_panel_weights <- unbalanced_panel_weights
    ### build initial pair structure
    cml_pair_type <- NULL
    if (!is.null(control_weights$cml_pair_type)) {
      cml_pair_type <- control_weights$cml_pair_type
    } else if (!is.null(Rprobit_obj$control$control_weights$cml_pair_type)) {
      cml_pair_type <- Rprobit_obj$control$control_weights$cml_pair_type
    }
    if (!is.null(cml_pair_type)) {
      Rprobit_obj$control$pairs_list <- build_cml_pairs(Rprobit_obj = Rprobit_obj, cml_pair_type = cml_pair_type)
    }


    ### check for unbalanced panel situation
    if (!is.null(control_weights$unbalanced_panel_weights)) {
      Rprobit_obj$control$control_weights$unbalanced_panel_weights <- control_weights$unbalanced_panel_weights

      use_panel_weights <- FALSE
      if (is.list(control_weights$unbalanced_panel_weights)) {
        use_panel_weights <- TRUE
      } else if (length(control_weights$unbalanced_panel_weights) == 1) {
        if (is.logical(control_weights$unbalanced_panel_weights)) {
          use_panel_weights <- control_weights$unbalanced_panel_weights
        }
      }

      #if (use_panel_weights) {
        panel_group_weights <- build_unbalanced_panel_weights(Rprobit_obj = Rprobit_obj)

        ### save estimated theta parameter, if applicable
        if (!is.null(attr(panel_group_weights, "theta"))) {
          Rprobit_obj$theta <- attr(panel_group_weights, "theta")
        }

        
        for (i_weight in 1:nrow(panel_group_weights)) {
          unbalanced_panel_weights[Rprobit_obj$data$Tp == panel_group_weights$n_obs[i_weight]] <- panel_group_weights$weights[i_weight]

        }
      #}
    }


    ### convert CML pair type to list
    if (is.null(cml_pair_type)) {
      cml_pair_type <- 0
    }
    if (is.numeric(cml_pair_type)) {
      cml_pair_type <- list(name = cml_pair_type)
    } else if (is.character(cml_pair_type)) {
      cml_pair_type <- list(name = cml_pair_type)
    }
    cml_pair_type$unbalanced_panel_weights <- unbalanced_panel_weights


    Rprobit_obj$control$control_weights$cml_pair_type <- cml_pair_type
    CML_pairs <- build_cml_pairs(Rprobit_obj = Rprobit_obj, cml_pair_type = cml_pair_type)

    Rprobit_obj$control$pairs_list <- CML_pairs

    out <- Rprobit_obj
  }

  return(out)
}
