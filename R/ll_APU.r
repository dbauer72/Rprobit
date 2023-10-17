#' implements the mixture of marginal and pairwise CML function
#' 
#' @description 
#' This function adds the marginal CML to the pairwise CML. 
#' 
#' @param theta
#' Parameter vector
#' @param data_obj
#' A \code{data_cl} object
#' @param control 
#' \code{control_cl} object controlling the evaluation of the likelihood. 
#' @param mod 
#' \code{\link{mod_cl}} object encoding the model.
#' 
#' @return 
#' A number containing the attributes "gradient" and "Hessian". 
#' 
#' @export

ll_APU <- function(theta, data_obj, mod, control) {
  out1 <- ll_macml_marginal(theta = theta, data_obj = data_obj, mod = mod, control = control)
  out2 <- ll_macml(theta = theta, data_obj = data_obj, mod = mod, control = control)
  
  gr1 =attr(out1,"gradient")
  H11 = attr(out1,"hessian1")
  H21 = attr(out1,"hessian2")
  
  gr2 =attr(out2,"gradient")
  H12 = attr(out2,"hessian1")
  H22 = attr(out2,"hessian2")
  
  out <- out1+out2
  attr(out,"gradient") <- gr1+gr2
  attr(out,"hessian1") <- H11+H12
  attr(out,"hessian2") <- H21+H22

  if (!is.null(attr(out1,'hessian'))){
    attr(out,"hessian") <- attr(out1,"hessian")+attr(out2,"hessian")
  }
  return(out)  
}
