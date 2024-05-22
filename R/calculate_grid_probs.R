#' Calculate the point masses at the grid points 
#'
#' @description
#' This function is a helper function used for the evaluation of the non-parametric mixing distributions.
#'
#' @param mod
#' A \code{mod_cl} object
#' @param theta
#' A parameter vector (size depends on mod_cl).
#'
#' @return
#' A real vector of probabilities
#'
#' @keywords internal

calculate_grid_probs <- function(mod,theta){
  if (inherits(mod,"mod_nonpara_splines_cl")){
    # spline calculations 
    K <- length(theta)
    
    # calculate spline basis
    spline_basis <- cal_spline_basis(mod$params,mod$knots)
    # regularization to avoid overly large or small values
    spline_out <- theta %*% spline_basis
    
    # calculate conditional distribution
    pi_est<- exp(spline_out)/(1+exp(spline_out))^2
    probs <- pi_est/sum(pi_est)
    
  } else {
    # full grid 
    th <- -6 + 12*exp(theta)/(1+exp(theta))
    
    # calculate conditional distribution
    pi_param <- c(0,th)
    probs<- exp(pi_param)/sum(exp(pi_param))
  }

  return (probs)  
}

