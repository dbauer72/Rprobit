#' Check if 'data_raw' is a valid 'data_raw'-object
#' @description Function that checks if 'data_raw' is a valid 'data_raw'-object.
#' @param data_raw \code{data_raw}-object
#' @param form \code{form}-object
#' @param ids vector specifying the columns in \code{data_raw} that identify a choice occasion uniquely
#' @param allvars \code{allvars}-object
#' @return Number of  alternatives found in \code{data_raw}
#' @keywords internal

check_data_raw <- function(data_raw, form, ids, allvars) {

  ### check if 'data_raw' is of the class 'data.frame'
  if (!is.data.frame(data_raw)) {
    stop("'data_raw' needs to be of the class 'data.frame'.")
  }

  ### check if 'data_raw' contains a named column with the choices
  choice <- all.vars(form)[1]
  if (!(choice %in% colnames(data_raw))) {
    stop("Dependent variable is not a column in 'data_raw'. Correct name?")
  }

  ### check if 'data_raw' contains all named columns that identify a choice occasion uniquely
  if (any(!ids %in% colnames(data_raw))) {
    stop("At least one identifier of 'ids' is not a column in 'data_raw'. Correct name(s)?")
  }

  ### compute number of alternatives in 'data_raw'
  alt <- check_allvars(data_raw = data_raw, allvars = allvars, alt = length(unique(data_raw[, choice])))

  return(alt)
}
