#' Read formula
#'
#' @description
#' Function that reads out different types of covariates from formula object 'form'.
#'
#' @param form
#' An object of class \code{formula} of the form
#' \code{choice ~ Type 1 | Type 2 | Type 3}, where:
#' - \code{choice} indicates the choice variable
#' - \code{Type 1} regressors are varying over alternatives and enter the model as is
#' - \code{Type 2} regressors are constant over alternatives and enter the model with alternative varying coefficients
#' (setting the coefficient for norming alternative to zero)
#' - \code{Type 3} regressors are alternative varying and enter the utilities with alternatives varying coefficients
#'
#' @return A list containing the following elements:
#' \item{allvars}{\code{allvars}-object (list of three entries listing the Type 1, 2 and 3 variables)}
#' \item{ASC}{Boolean, indicating whether ASCs are estimated or not}
#'
#' @keywords internal

read_formula <- function(form){

  varraw <- trimws(strsplit(as.character(form)[3], split="|", fixed = TRUE)[[1]])
  allvars <- lapply(strsplit(varraw, split="+" , fixed=TRUE), trimws)
  var_levels <- length(allvars)
  ASC <- TRUE

  if(var_levels==1){
    allvars[[2]] = 0
    allvars[[3]] = 0
  }

  if(var_levels==2){
    allvars[[3]] = 0
  }

  if(var_levels > 1){

    ### check if there is a 0 AND 1 in the second spot
    if (any(allvars[[2]] %in% 0) && any(allvars[[2]] %in% 1)) {
      stop("A '0' and a '1' were included in the second spot of the formula, which leads to a contradiction. Please chose either '0' or '1' to decide if ASCs should be suppressed or not.")
    }

    ### check if ASCs should be suppressed
    if (any(allvars[[2]] %in% 0)) {
      ASC <- FALSE
    }

    ### check if ASCs should be included
    if (any(allvars[[2]] %in% 1)) {
      ASC <- TRUE
    }
  }

  ### remove redundant 0's and 1's in any spot
  for (i in 1:3) {
    if (any(allvars[[i]] %in% c(0, 1))) {
      allvars[[i]] <- allvars[[i]][!(allvars[[i]] %in% c(0, 1))]
      if (length(allvars[[i]]) == 0) {
        allvars[[i]] <- 0
      }
    }
  }

  return(list("allvars" = allvars, "ASC" = ASC))
}
