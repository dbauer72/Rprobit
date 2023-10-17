#' Find the empirical likelihood value at CMLE for model selection 
#' 
#' @description 
#' This function estimates the Lagrangian parameter in the empirical likelihood framework and returns the corresponding EL value.
#' 
#' @param Rprobit_obj 
#' An \code{\link{Rprobit_cl}} object.
#' @return 
#' A real number, the maximum of the empirical likelihood criterion function. 
#' 
#' @export
calc_emp_like <- function(Rprobit_obj) {
  
  # clone to avoid sideeffects. 
  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  
  # set control to the evaluation of scores needed for EL. 
  Rprobit_o$control$el = 1 
  data_tr <- Rprobit_o$data$clone()
  data_tr$data <- substract_choice_regressor_from_data(Rprobit_o$data$data)
  
  # run  CML function to calculate score and derivative
  out <- ll_macml(theta = Rprobit_o$theta, data_obj = data_tr, mod = Rprobit_o$mod, control = Rprobit_o$control)
  
  # retrieve score function 
  score = attr(out,"gradEL")
  # grad_score = attr(out,"HessEL")
  
  # convert to matrix 
  N = length(score)
  P = length(score[[1]])
  sn = matrix(0,P,N)
  for (n in 1:N){
    sn[,n] = score[[n]]
  }

  logp <- function(x,N){
    if (x>1/N){
      val = log(x)
    } else {
      slope = N
      val = -log(N)+slope*(x-1/N)-N^2*(x-1/N)^2
    }
    return(val)
  }
  
  # set up EL and its gradient 
  eval_lE <- function (psi, sn){
    psi_sn <- as.vector(psi %*% sn)
    N = dim(sn)[2]
    #psi_sn <- psi_sn /(1+exp(-N*min(psi_sn)))
    #grad = rep(0,length=length(psi))
    if (min(psi_sn)<=(1/N-1)){
      lE = 100000000
    #  grad = -psi
    } else {
      #lE <- 2*sum(unlist(lapply(1+ psi_sn,logp,N=N)))
      #lE <- 2*sum(log(1+psi_sn))
      gn <- sn %*% (1/(1+psi_sn))
      lE <- sum ( gn^2)              
    #  grad = -2*(lapply(psi_sn,sn %*% (1/(1+psi_sn)))
    }

    #attr(lE,"gradient") = grad;
    return (lE)
  }
  

  
  psi_init = rep(0,P)
  # perform numerical optimization

  out <- stats::nlm(f                  = eval_lE,
                    p                  = psi_init,
                    sn                 = sn,
                    print.level        = 0)
  
  
  # retrieve optimal value. 
  psi_est <- out$estimate
  psi_sn <- as.vector(psi_est %*% sn)
  psi_sn <- pmax(psi_sn,rep(-1+0.0000000001,length(psi_sn)))
  el_val <- 2*sum(log(1+psi_sn))
  
  return (el_val)
}
