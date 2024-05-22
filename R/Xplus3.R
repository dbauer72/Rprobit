#' evaluate (x)_+^3 for splines 
#'
#' @description
#' This function is a helper function used for the evaluation of the spline calculations.
#'
#' @param x
#' A real value
#'
#' @return
#' A real value.
#'
#' @keywords internal

Xplus3 <- function(x){
  
  xp3 <- x*0
  for (j in 1:length(x)){
    xp3[j] = (max(x[j],0))^3
  }
  return (xp3)
}

