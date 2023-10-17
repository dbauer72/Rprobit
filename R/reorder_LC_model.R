#' Build a covariance matrix and calculates the standard deviations
#' @description Function that builds a covariance matrix and calculates the standard deviations using the Delta method.
#' @param Rprobit_obj
#' A \code{\link{Rprobit_cl}} object.
#' @return The reordered \code{\link{Rprobit_cl}} object.
#' @keywords internal

reorder_LC_model <- function(Rprobit_obj) {
  object <- Rprobit_obj$clone(deep = TRUE) # clone object to avoid side effects.
  mod <- object$mod
  theta <- object$theta
  
  # check, if more than one class -> if not, do nothing. 
  if (!is.null(mod$num_class)){
    # extract number of classes 
    num_class <- mod$num_class
    param_one <- mod$lthb + mod$lthO + mod$lthL

    # construct matrix of params and indices
    theta_mat <- matrix(theta[c(1:(param_one*num_class))], nrow = param_one)
    ind_mat <- matrix(c(1:(param_one*num_class)),nrow = param_one)
    
    # calculate pi mixing probs 
    pi_eta = c(0,theta[num_class*param_one + c(1:(num_class-1))])
    pi = exp(pi_eta)/sum(exp(pi_eta))
    
    # sort probs descending 
    ind_sort <- order(pi, decreasing = TRUE)

    # reorder matrix of parameters 
    theta_par <- matrix(theta_mat[,ind_sort],ncol=1)
    ind_vec <- matrix(ind_mat[,ind_sort],ncol=1)
    
    theta[c(1:(param_one*num_class))] <- theta_par
    
    # reorder gamma parameters for mixing distribution 
    gamma <- pi_eta[ind_sort]
    gamma <- gamma - gamma[1]
    
    theta[num_class*param_one + c(1:(num_class-1))] <- gamma[-1]

    object$theta <- theta 
    
    # adjust variance and Hessian. 
    if (any(is.na(object$vv)==FALSE)){
      # parameters for classes 
      ind_old = c(1:(param_one*num_class))
      object$H[ind_old,ind_old] <- object$H[ind_vec,ind_vec]
      object$J[ind_old,ind_old] <- object$J[ind_vec,ind_vec]
      object$vv[ind_old,ind_old] <- object$vv[ind_vec,ind_vec]
      if (any(is.na(object$vv2)==FALSE)){
        object$vv2[ind_old,ind_old] <- object$vv2[ind_vec,ind_vec]
      }
      
      # parameters for mixing probabilities 
      Trafo <- matrix(0,num_class-1,num_class-1)
      P <- diag(num_class)
      P <- P[ind_sort,]
      P <- P[,-1]
      T2 <- diag(num_class)
      T2 <- T2 - matrix(1,num_class,1) %*% T2[,1]
      Trafo <- T2 %*% P
      Trafo <- Trafo[-1,]
      tTrafo <- t(Trafo)      
      ind_gam <- num_class*param_one + c(1:(num_class-1))
      # write results in matrices
      object$H[ind_old,ind_gam] <- object$H[ind_old,ind_gam] %*% tTrafo 
      object$H[ind_gam,ind_old] <- Trafo %*% object$H[ind_gam,ind_old]
      object$H[ind_gam,ind_gam] <- Trafo %*% object$H[ind_gam,ind_gam] %*% tTrafo 

      object$J[ind_old,ind_gam] <- object$J[ind_old,ind_gam] %*% tTrafo 
      object$J[ind_gam,ind_old] <- Trafo %*% object$J[ind_gam,ind_old]
      object$J[ind_gam,ind_gam] <- Trafo %*% object$J[ind_gam,ind_gam] %*% tTrafo 

      object$vv[ind_old,ind_gam] <- object$vv[ind_old,ind_gam] %*% tTrafo 
      object$vv[ind_gam,ind_old] <- Trafo %*% object$vv[ind_gam,ind_old]
      object$vv[ind_gam,ind_gam] <- Trafo %*% object$vv[ind_gam,ind_gam] %*% tTrafo 
      
      if (any(is.na(object$vv2)==FALSE)){
        object$vv2[ind_old,ind_gam] <- object$vv2[ind_old,ind_gam] %*% tTrafo 
        object$vv2[ind_gam,ind_old] <- Trafo %*% object$vv2[ind_gam,ind_old]
        object$vv2[ind_gam,ind_gam] <- Trafo %*% object$vv2[ind_gam,ind_gam] %*% tTrafo       }
    }
    
  }  
  
  return(object)
}
