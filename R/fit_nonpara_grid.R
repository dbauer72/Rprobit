#' Fit a Probit Model using non-parametric mixing distribution based on a fixed grid.
#' The parameter for the first regressor is mixed non-parametrically, the remaining parameters are fixed and contained in 
#' the grid elements.  
#' 
#' @description 
#' This function fits a mixed probit model using a non-parametric mixing with fixed grid and estimated mixing 
#' probabilities.  
#' DOES NOT USE PARALLELIZATION!!!
#' 
#' @param data_tr 
#' data set already preprocessed
#' @param mod 
#' \code{\link{mod_cl}} object
#' @param control
#' \code{\link{control_cl}} object
#' @param cml_pair_type 
#' An integer. By default, \code{cml_pair_type = 1}. -1 triggers independence likelihood, 0 for FP, 1 for AP. 
#' 
#' 
#' @return 
#' A \code{\link{Rprobit_cl}} object.
#' 
#' @export
#' 

fit_nonpara_grid <- function(data_tr, mod, control, cml_pair_type  = 1) {

  probs <- choice_probs_nonpara(data_tr, mod, control, cml_pair_type)
  #grid_points <- mod$params 
  
  # initialize 
  prs <- apply(probs,2,mean)
  #prs <- prs - prs[1]
  pi <- prs/sum(prs)
  theta <- pi
  
  # evaluate negative log-likelihood of non-parametric model  
  eval_ll_grad <- function(pi,probs){
    tol = 0.000000001
    K <- dim(probs)[2]
    N <- dim(probs)[1]
  
    # criterion function 
    pri = probs %*% pi
    for (j in 1:N){
      probs[j,] <- -probs[j,]/max(pri[j],tol)
    }
    grad <- apply(probs,2,sum)
  
    return (grad)
  }

  # evaluate negative log-likelihood of non-parametric model  
  eval_ll <- function(pi,probs){
    tol = 0.000000001
    K <- dim(probs)[2]
    N <- dim(probs)[1]
  
    # criterion function 
    pri = probs %*% pi
    LN <- (-1)*sum( log(pri))
  
    # gradient 
    for (j in 1:N){
      probs[j,] <- -probs[j,]/max(pri[j],tol)
    }
    grad <- apply(probs,2,sum)
  
    return( list( "objective" = LN,"gradient"  = grad ))
  }
  
  eval_g_eq <- function( x ) {
    constr <- sum(x)-1
    grad <- matrix(1,1,length(x))
    
    return( list( "constraints"=constr, "jacobian"=grad ) )
  }
  
  eval_ll_x <- function(x){ eval_ll(x,probs=probs)}
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-7 )  
  opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", "xtol_rel" = 1.0e-7,  "maxeval" = 10000, "local_opts" = local_opts )
  
  out_opt <- nloptr::nloptr(x0=theta,eval_f=eval_ll_x,lb = rep(0,length(theta)),eval_g_eq=eval_g_eq, opts = opts)

  # provide return values 
  pi <- out_opt$solution 
  sprob <- probs %*%  out_opt$solution
  return (list(lprob = sprob, pi = pi))
} 
