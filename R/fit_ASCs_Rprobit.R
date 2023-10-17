#' Fit ASCs for an estimated Probit Model
#' 
#' @description 
#' This function adjusts the ASCs for an estimated probit model such that 
#' mean predicted frequencies and observed frequencies coincide.
#' 
#' @param Rprobit_obj 
#' An \code{\link{Rprobit_cl}} object.
#'
#' @return 
#' A \code{\link{Rprobit_cl}} object.
#' 
#' @export

fit_ASCs_Rprobit <- function(Rprobit_obj) {
  
  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  mod_orig <- Rprobit_o$mod$clone(deep=TRUE)
  data_raw_orig = Rprobit_o$data_raw$clone(deep=TRUE)
  
  # find variables for ASCs 
  vars <- Rprobit_o$vars 
  
  ind <- which(startsWith(vars,"ASC"))
  # find parameters for ASCs 
  Hb <- mod_orig$Hb
  
  param_ind <- ind*0
  for (j in 1:length(ind)){
    param_ind[j] <- which(Hb[ind[j],]==1) # This assumes that ASCs enter withouth any linear restriction. 
  }
  
  theta_init <- Rprobit_o$theta[param_ind]
  theta <- Rprobit_o$theta
  theta[param_ind] <- 0 
  b <- Hb %*% theta[1:dim(Hb)[2]]
  
  theta_O <- theta[(mod_orig$lthb+1):(mod_orig$lthb+mod_orig$lthO)]
  theta_sig <- theta[(mod_orig$lthb+mod_orig$lthO+1):(mod_orig$lthb+mod_orig$lthO+mod_orig$lthL)]
  
  # set up new model structure
  mod <- mod_orig$clone(deep = TRUE)
  mod$fb <- b 
  mod$Hb <- Hb[,ind]
  
  if (mod$lthO>0){
    mod$fO <- mod$fO + mod$HO %*% theta_O
    mod$HO <- matrix(0,dim(mod$fO)[1],0)
    mod$lRE <- 0
  }

  if (mod$lthL>0){
    mod$fL <- mod$fL + mod$HL %*% theta_sig
    mod$HL <- matrix(0,dim(mod$fL)[1],0)
  }
  
  Rprobit_o$mod <- mod
  Rprobit_o$theta_0 <- theta_init
  Rprobit_o$theta <- theta_init
  
  # run estimation with new model structure
  cal_diff_pred <- function(theta,Rprobit_o){
    Rprobit_o$theta <- theta
    pro <- predict_Rprobit(Rprobit_o)
    err <- sum(abs(diff(pro)))
    return (err)
  }
  
  out <- stats::nlm(f                  = cal_diff_pred,
                    p                  = theta_init,
                    Rprobit_o          = Rprobit_o)
  
  #Rprobit_o <- fit_Rprobit(Rprobit_o, init_method = "theta", cml_pair_type=0)
  
  # fill in new parameters for old model structure 
  theta[param_ind] <- out$estimate #Rprobit_o$theta[1:length(param_ind)]
  Rprobit_o$theta <- theta
  Rprobit_o$mod <- mod_orig 
  
  # re-evaluate likelihood and variance 
  
  time_1                    <- Sys.time()
  #if(Rprobit_o$control$approx_method!="TVBS"){

  if(is.null(Rprobit_o[["data"]]) & !is.null(Rprobit_o[["data_raw"]])){
    data_from_data_raw_model <- function(Rprobit_o){
      read_formula_out  <- read_formula(Rprobit_o$form)
      ASC               <- read_formula_out$ASC
      allvars           <- read_formula_out$allvars
      choice            <- all.vars(Rprobit_o$form)[1]
      norm_alt          <- Rprobit_o$info$setup_input$norm_alt
      alt               <- Rprobit_o$mod$alt
      
      data_raw_to_data_out = data_raw_to_data(data_raw = Rprobit_o$data_raw, allvars = allvars, choice = choice, re = Rprobit_o$re, norm_alt = norm_alt, alt = alt, ASC = ASC)

      return(data_raw_to_data_out)
    }
    
    Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
  }
  data_tr <- Rprobit_o$data$clone()
  data_tr$data <- substract_choice_regressor_from_data(Rprobit_o$data$data)
  
  neg_ll_fit                <- ll_macml(theta = Rprobit_o$theta, data_obj = data_tr, mod = mod_orig, control = Rprobit_o$control)
  
  time_2                    <- Sys.time()
  Rprobit_o$ll    <- (-1)*(neg_ll_fit)
  Rprobit_o$grad <- attr(neg_ll_fit,"gradient")

  
  if (Rprobit_o$control$hess == TRUE){
    Rprobit_o$H               <- attr(neg_ll_fit, "hessian")
  } else {
    Rprobit_o$H <- attr(neg_ll_fit, "hessian2")
  }
  if (Rprobit_o$control$probit == TRUE){
    Rprobit_o$J               <- attr(neg_ll_fit,"hessian2")
  } else {
    Rprobit_o$J               <- attr(neg_ll_fit,"hessian1")
  }

  ### save individual contributions to gradient and hessian
  if ("gradEL" %in% names(attributes(neg_ll_fit))) {
    Rprobit_o$gradEL <- attr(neg_ll_fit, "gradEL")
  }
  if ("HessEL" %in% names(attributes(neg_ll_fit))) {
    Rprobit_o$HessEL <- attr(neg_ll_fit, "HessEL")
  }

  ### if covariance matrix H is not pos.definite, adjust
  H_fixed                   <- adjust_variance_matrix(Rprobit_o$H)
  if (H_fixed$ok == FALSE){
    message("Calculated variance matrix H needed to be adjusted due to negative eigenvalue(s)!")
  }
  #inverted_H                <- solve(H_fixed$mat)
  inverted_H                <- H_fixed$mat_inv

  ### if covariance matrix J is not pos.definite, adjust
  S = adjust_variance_matrix(Rprobit_o$J)
  if (S$ok == FALSE){
    message("Calculated variance matrix J needed to be adjusted due to eigenvalue(s) close to zero!")
  }
  if (Rprobit_o$control$probit == TRUE) {
    # Rprobit_o$vv              <- solve(S$mat)
    Rprobit_o$vv <- S$mat_inv
  } else {
    Rprobit_o$vv <- inverted_H %*% S$mat %*% inverted_H
  }
  time_2 <- Sys.time()

  ### save time required for hessian calculation in Rprobit_o$info$hess_time
  if (!is.null(Rprobit_o$info)) {
    Rprobit_o$info$hess_time <- difftime(time_2, time_1, units = "secs")
  } else {
    Rprobit_o$info <- list(hess_time = difftime(time_2, time_1, units = "secs"))
  }


  ### stop using multiple clusters
  # stop_parallel(Rprobit_o = Rprobit_o, use_parallel = use_parallel)

  ### save Rprobit Package version in Rprobit_o$info$Rprobit_version
  R_version_nr <- getRversion()
  Rprobit_version_nr <- utils::packageVersion(pkg = "Rprobit")
  if (!is.null(Rprobit_o$info)) {
    Rprobit_o$info$Rprobit_version <- Rprobit_version_nr
    Rprobit_o$info$R_version <- R_version_nr
  } else {
    Rprobit_o$info <- list(Rprobit_version = Rprobit_version_nr,
                           R_version = R_version_nr)
  }
  rm(Rprobit_version_nr, R_version_nr)


  ### Remove list object data from Rprobit_o to reduce storage space
  Rprobit_o$data <- NULL
  Rprobit_o$mod <- mod_orig
  Rprobit_o$data_raw <- data_raw_orig

  return(Rprobit_o)
}
