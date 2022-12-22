#' Build Info Part of the Model
#'
#' @description
#' This function builds the info part of the model.
#'
#' @param Rprobit_obj
#' An \code{\link{Rprobit_cl}} object.
#' @param data_raw
#' \code{data_raw}-object
#' @inheritParams read_formula
#' @param name
#' string, Name for the model
#' @param ids
#' string, vector of variable names of the data object, identifying the variables that identify a decision maker
#' @param description
#' string, Description for the model
#' @param control
#' List, contains control variables for estimation procedure
#' @param info
#' List, contains information on the model and estimation performance
#'
#' @return
#' An \code{\link{Rprobit_cl}} object.
#'
#' @keywords internal

build_info <- function(
    Rprobit_obj, data_raw, form, ids, name, description, control, info = NULL
  ) {

  ### initialize info list, if not existent
  if (is.null(info)) {
    info <- list()
  }

  ### save the name of the columns which identify the individuals
  info$ids <- ids

  ### build the name of the model or retrieve it
  if (is.null(name)) {
    if (is.null(info$name)) {
      name <- "Discrete Choice Model"
      if (!is.null(control$probit)) {
        if (control$probit) {
          name <- paste("Probit", name)
        } else {
          name <- paste("MACML", name)
        }
      }
      if (is.null(data_raw)) {
        name <- paste("Simulated", name)
      }
    } else {
      name <- info$name
    }
  }

  ### build the drscription of the model or retrieve it
  if (is.null(description)) {
    if (is.null(info$description)) {
      description <- "Discrete Choice Model"

      Tp_summary <- c(max(Rprobit_obj$data_raw$Tp), min(Rprobit_obj$data_raw$Tp))
      if (Tp_summary[1] == Tp_summary[2]) {
        description_tp <- paste("T=", Tp_summary[1], " obsevations each", sep = "")
      } else {
        description_tp <- paste("between T=", Tp_summary[2], " and T=", Tp_summary[1], " obsevations each", sep = "")
      }
      description_add <- paste("with N=", Rprobit_obj$data_raw$N, " individuals with ", description_tp, " and ", Rprobit_obj$mod$alt, " alternatives", sep = "")

      description <- paste(description, " ", description_add, sep = "")


      if (!is.null(control$probit)) {
        if (control$probit == TRUE) {
          description <- paste("Probit", description)
        } else {
          description <- paste("MACML", description)
        }
      }

      if (!is.null(control$approx_method)) {
        if (control$approx_method == "SJ") {
          description <- paste(description, "using an SJ-approximation for the Gaussian distributions")
        }
        if (control$approx_method == "TVBS") {
          description <- paste(description, "using a TVBS-approximation for the Gaussian distributions")
        }
      }

      if (is.null(data_raw)) {
        description <- paste("Simulated ", description, ". True model parameters are theta_0=(", paste(Rprobit_obj$theta, collapse = ", "), ").", sep = "")
      } else if (!is.null(form)) {
        description_add <- paste(form)
        description <- paste(description, "using the model formula ", description_add[2], description_add[1], description_add[3])
      }
    } else {
      description <- info$description
    }
  }


  ### add the name and description to the info
  info$name <- name
  info$description <- description


  return(info)
}
