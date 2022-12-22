#' Check alignment between 'allvars' and 'data_raw'
#' @description Function that checks the alignment between 'allvars' and 'data_raw' and returns the found number of alternatives.
#' @param data_raw \code{data_raw}-object
#' @param allvars \code{allvars}-object
#' @param alt Number of choice alternatives
#' @return Number of choice alternatives
#' @keywords internal

check_allvars <- function(data_raw, allvars, alt) {

  ### initial number of alternatives
  alt0 <- alt

  ### check whether alternative specific covariates (first or third type) are present in 'data_raw' for all alternatives
  for (i_var in c(allvars[[1]], allvars[[3]])) {
    if (i_var != 0) {
      var_occ <- length(grep(paste0("^", i_var, "_"), colnames(data_raw)))
      if (var_occ < alt) {
        stop(paste("The following alternative specific covariate is not supplied for all alternatives:", i_var))
      }
      if (var_occ > alt) {
        warning(paste("The number of alternatives is increased by", alt - var_occ, "because the following alternative specific covariate is supplied for additional alternatives:", i_var))
        alt <- var_occ
      }
    }
  }

  ### check whether covariates of second type are present in 'data_raw'
  for (i_var in allvars[[2]]) {
    if (i_var != 0) {
      var_occ <- length(grep(paste0("^", i_var, "$"), colnames(data_raw)))
      if (var_occ == 0) {
        stop(paste("The following covariate is missing from the data set:", i_var))
      }
    }
  }

  ### repeat 'check_allvars' if additional, unchosen alternatives have been discovered from the alternative specific variables
  if (alt > alt0) {
    check_allvars(data_raw, allvars, alt)
  } else {
    return(alt)
  }
}
