#' Initialize fit_LC_Rprobit
#'
#' @description
#' Function that initializes the estimation routine for latent classes model.
#'
#' @param Rprobit_obj
#' an \code{\link{Rprobit_cl}} object. 
#' @param mod
#' An \code{\link{mod_cl}} object.
#' @param control
#' An \code{\link{mod_cl}} object. 
#'
#' @return
#' a vector of model parameters
#'
#' @export

init_LC_Rprobit <- function(Rprobit_obj,mod,control) {

  K <- mod$num_class
  # estimate joint model without latent classes 
  Rp <- fit_Rprobit(Rprobit_obj, init_method = "theta")
  Rprobit_o = Rprobit_obj$clone(deep = TRUE)
  data_tr = Rprobit_o$data$clone()
  data_tr$data <- substract_choice_regressor_from_data(Rprobit_o$data$data)
  
  out <- ll_macml(theta = Rp$theta, data_obj =data_tr, mod = Rprobit_obj$mod, control = Rprobit_obj$control)

  param_one <- length(Rp$theta)
  
  # collect gradient data into matrix 
  grad <- t(matrix(unlist(attr(out,"gradEL")),nrow = length(Rp$theta)))

  # grad is a wide matrix: N x npar x: standardize
  # due to optimization row means are zero. 
  N = dim(grad)[1]
  V = t(grad) %*% grad/N
  L = chol(V)
  iL = solve(L)
  
  # normalize gradient matrix. 
  grad <- grad %*% iL 

  # cluster gradient matrix  
  km <- stats::kmeans(grad,centers = K,nstart = 5)
  
  # for each cluster re-estimate the model 
  theta = c()
  for (k in 1:K){
    ind = which(km$cluster ==k)
    Rp_comp <- Rp$clone(deep = TRUE)
    data_sub <- Rprobit_obj$data$clone(deep = TRUE)
    data_sub$set_data(Rprobit_obj$data$data[ind])
    Rp_comp$data = data_sub
    Rp_comp$data_raw$set_df(Rp_comp$data_raw$df[c(1:(5*length(ind))),])
    Rp_comp <- fit_Rprobit(Rp_comp, init_method = "theta")
    theta <- c(theta,Rp_comp$theta)
  }
  
  # add parameters for pi_est. 
  freq <- table(km$cluster)
  pi <- freq/sum(freq)
  gamma_est <- log(pi[-1]/pi[1])
  theta <- c(theta,gamma_est)
  
  return (theta) 
}
 
