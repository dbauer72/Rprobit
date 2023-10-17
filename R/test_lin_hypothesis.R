#' Test linear hypotheses on the parameter vector using empirical likelihood and adjusted likelihood ratios. 
#' 
#' @description 
#' This function uses the estimate without restrictions to obtain an estimate under the restricted CML.
#' Subsequently different tests for a linear hypothesis on the parameter vector is performed 
#' using the empirical likelihood ratio as well as adjustments to the CML ratio. 
#' 
#' @param Rprobit_obj 
#' A \code{\link{Rprobit_cl}} object.
#' @param hypothesis
#' Either a \code{\link{hypothesis_cl}} object or a vector for H_0: theta = theta_0. 
#'
#' @return 
#' A real matrix where each row has two entries corresponding to the unrestricted and the restricted estimate of
#' + CML ratio
#' + empirical likelihood ratio
#' + AIC
#' + BIC
#' 
#' @export
test_lin_hypothesis <- function(Rprobit_obj,hypothesis) {
  
  # initialize criterion matrix 
  crits = matrix(0,8,3)
  
  # obtain characteristics of estimation
  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  data <- Rprobit_obj$data$clone()
  
  # calculate emp  
  emp_un <- calc_emp_like(Rprobit_o)
  
  crits[1,1] = emp_un
  crits[2,1] = Rprobit_o$ll
  crits[3,1] = Rprobit_o$AIC()
  crits[4,1] = Rprobit_o$BIC()
  H = Rprobit_o$H
  J = Rprobit_o$J
  dfs = sum(diag( solve(H) %*% J))
  crits[5,1] = -2*Rprobit_o$ll + 2*dfs
  crits[6,1] = -2*Rprobit_o$ll + log(Rprobit_o$data_raw$N)*dfs
  
  mod_orig <- Rprobit_o$mod$clone(deep=TRUE)
  
  H0_theta0 = 0; 
  
  if (inherits(hypothesis,"hypothesis_cl")==TRUE){

    newparam = dim(hypothesis$Hb)[2]+ dim(hypothesis$HO)[2]+dim(hypothesis$HL)[2]
    if (newparam>0){ # parameters left under H_0.
      Rprobit_r <- Rprobit_o$clone(deep = TRUE) # new object for optimization of restricted model. 
      thetahat <- Rprobit_o$theta
      if (mod_orig$lthb>0){
        thetahat_beta <- thetahat[1:mod_orig$lthb]
        mod_orig$fb = mod_orig$fb + mod_orig$Hb %*% hypothesis$fb
        mod_orig$Hb = mod_orig$Hb %*%hypothesis$Hb
      } else {
        thetahat_beta = NULL
      }
      if (mod_orig$lthO>0){
        thetahat_omega <- thetahat[mod_orig$lthb+(1:mod_orig$lthO)]
        mod_orig$fO = mod_orig$fO + mod_orig$HO %*% hypothesis$fO
        mod_orig$HO = mod_orig$HO %*%hypothesis$HO
      } else {
        thetahat_omega = NULL
      }
      if (mod_orig$lthL>0){
        thetahat_sigma <- thetahat[mod_orig$lthb+mod_orig$lthO + c(1:mod_orig$lthL)]
        mod_orig$fL = mod_orig$fL + mod_orig$HL %*% hypothesis$fL
        mod_orig$HL = mod_orig$HL %*% hypothesis$HL  
      } else {
        thetahat_sigma = NULL
      }
      # adjust object to incorporate hypothesis
      
      Rprobit_r$mod <- mod_orig
    
      # define function for pseudo-inverse
      pinv <- function(R){
        Rdag= solve( t(R) %*% R) %*% t(R)
        return(Rdag)
      }
      
      # obtain estimate 
      if (dim(hypothesis$Hb)[2]>0){
        tilde_beta = pinv(hypothesis$Hb) %*% (thetahat_beta- hypothesis$fb)
      } else { 
        tilde_beta = NULL 
      }
      if (dim(hypothesis$HO)[2]>0){
        tilde_omega = pinv(hypothesis$HO) %*% (thetahat_omega- hypothesis$fO)
      } else {
        tilde_omega = NULL
      }
      if (dim(hypothesis$HL)[2]>0){
        tilde_sigma = pinv(hypothesis$HL) %*% (thetahat_sigma- hypothesis$fL)
      } else {
        tilde_sigma = NULL
      }
    
      tildetheta = c(tilde_beta,tilde_omega,tilde_sigma)
      # restrictions matrix
      R = matrix(0,length(thetahat),length(tildetheta))
      if (!is.null(tilde_beta)){
        R[1:length(thetahat_beta),1:length(tilde_beta)] = hypothesis$Hb
      }
      if (!is.null(tilde_omega)){
        R[length(thetahat_beta)+(1:length(thetahat_omega)),length(tilde_beta)+(1:length(tilde_omega))] = hypothesis$HO
      }
      if (!is.null(tilde_sigma)){
        R[length(thetahat_beta)+length(thetahat_omega)+(1:length(thetahat_sigma)),length(tilde_beta)+length(tilde_omega)+(1:length(tilde_sigma))] = hypothesis$HL
      }
    
      Rprobit_r$theta = tildetheta
      Rprobit_r$theta_0 = tildetheta
      Rprobit_r$data <- data
      Rprobit_r$control$control_nlm$typsize <- rep(1,length(tildetheta))
    
      Rprobit_r <- fit_Rprobit(Rprobit_r)
      
      # derive different criterion functions 
      tildetheta = Rprobit_r$theta
      
      # vector of fixed entries of theta
      f_all = as.matrix(c(hypothesis$fb,that_omega = hypothesis$fO,that_sigma = hypothesis$fL),ncol=1)
      
#      # obtain estimate 
#      # dimensions 
#      pb = dim(hypothesis$Hb)[2]
#      pO = dim(hypothesis$HO)[2]
#      pL = dim(hypothesis$HL)[2]
#      
#
#      
#      if (pb>0){
#        that_beta = that_beta + hypothesis$Hb %*% tildetheta[1:pb]
#      } 
#      if (pO>0){
#        that_omega = that_omega + hypothesis$HO %*% tildetheta[pb+(1:pO)]
#      }
#      if (pL>0){
#        that_sigma = that_sigma + hypothesis$HL %*% tildetheta[pb+pO+(1:pL)]
#      }
#      
#      Rprobit_o$theta = c(that_beta,that_omega,that_sigma)
      Rprobit_o$data <- data
      
      crit_emp_like <- function(tth,Rprobit_o,f_all,R_all){
        Rprobit_o$theta = f_all + R_all %*% tth
        fval <- calc_emp_like(Rprobit_o)

        return(fval)
      } 
        
      # optimize the function. 
      out <- stats::nlm(f                  = crit_emp_like,
                        p                  = tildetheta,
                        Rprobit_o          = Rprobit_o,
                        f_all              = f_all,
                        R_all              = R,
                        print.level        = 2)
      
      emp_r <- out$minimum
      
      # adjustment for CML ratios.
      Rdag = pinv(R)
      HmRHR = solve(H)- t(Rdag) %*% solve( Rdag %*% H %*% t(Rdag)) %*% Rdag 
      lambda = Re(eigen( HmRHR %*% J)$values)
      p = abs(diff(dim(R)))
      
    } else {  
      H0_theta0 = 1; # H0: theta = theta_0. No params left. 
      
      Rprobit_r <- Rprobit_o$clone(deep = TRUE) # new object for optimization of restricted model. 
      theta = c(hypothesis$fb,hypothesis$fO,hypothesis$fL)
      Rprobit_r$theta = theta
      
      Rprobit_r$control$el = 1 
      data_tr <- Rprobit_r$data$clone()
      data_tr$data <- substract_choice_regressor_from_data(Rprobit_r$data$data)
      
      # run  CML function to calculate score and derivative
      out <- ll_macml(theta = Rprobit_r$theta, data_obj = data_tr, mod = Rprobit_r$mod, control = Rprobit_r$control)
      Rprobit_r$ll <- -out
      
      
      
      # calculate emp  
      emp_r <- calc_emp_like(Rprobit_r)
      
      # adjustment factors 
      lambda = Re(eigen( solve(Rprobit_r$H) %*% Rprobit_r$J)$values)
      p = length(theta)
    }
    


  } else { # H_0: theta = theta_0. 
    # obtain characteristics of estimation
    H0_theta0 = 1; # H0: theta = theta_0. No params left. 
    
    Rprobit_r <- Rprobit_o$clone(deep = TRUE) # new object for optimization of restricted model. 
    Rprobit_r$theta = hypothesis
  
    Rprobit_r$control$el = 1 
    data_tr <- Rprobit_r$data$clone()
    data_tr$data <- substract_choice_regressor_from_data(Rprobit_r$data$data)
    
    # run  CML function to calculate score and derivative
    out <- ll_macml(theta = Rprobit_r$theta, data_obj = data_tr, mod = Rprobit_r$mod, control = Rprobit_r$control)
    Rprobit_r$ll <- -out
 
    # calculate emp  
    emp_r <- calc_emp_like(Rprobit_r)
    
    # adjustment factors 
    lambda = Re(eigen( solve(Rprobit_r$H) %*% Rprobit_r$J)$values)
    p = length(hypothesis)

  }

  # fill in results 
  crits[1,2] = emp_r
  crits[2,2] = Rprobit_r$ll
  crits[3,2] = Rprobit_r$AIC()
  if (H0_theta0 ==1){
    crits[3:6,2] = -2*Rprobit_r$ll
  } else {
    crits[3,2] = Rprobit_r$AIC()
    crits[4,2] = Rprobit_r$BIC()
    dfs = sum(diag( solve(Rprobit_r$H) %*% Rprobit_r$J))
    crits[5,2] = -2*Rprobit_r$ll + 2*dfs
    crits[6,2] = -2*Rprobit_r$ll + log(Rprobit_r$data_raw$N)*dfs
  }
  

  
  # last column is differences 
  crits[1,3] = crits[1,2]-crits[1,1]
  crits[2:6,3] = crits[2:6,1]-crits[2:6,2]
  
  # 7_ CML1 weights by average eigenvalue. 
  crits[7,] = crits[2,]
  crits[7,3] = crits[7,3]/sum(lambda)*p

  # 8_ CML2 weights by average eigenvalue. 
  crits[8,] = crits[2,]
  crits[8,3] = crits[8,3]/sum(lambda^2)*sum(lambda)
  crits[8,1] = sum(lambda)^2/sum(lambda^2)
  
  rownames(crits) <- c("EL","L_CML","AIC","BIC","CAIC","CBIC","C1LCML","C2LCML")
  colnames(crits) <- c("unr.","restr.","Diff")
  
  # return matrix
  return (crits)
}
