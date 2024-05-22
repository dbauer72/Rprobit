#' Create Rprobit-object from legacy mod_macml object
#'
#' @description
#' This function creates an \code{\link{Rprobit_cl}} object from a legacy mod_macml object from previous package iterations.
#'
#' @param mod_object
#' \code{mod_macml}-object
#'
#' @return
#' An \code{\link{Rprobit_cl}} object.
#'
#' @keywords internal

Rprobit_cl_from_modMACML <- function(mod_object){
  out <- Rprobit_cl$new(data_raw =  mod_object$data_raw, 
                                  data =      mod_object$dara, 
                                  form =      mod_object$form, 
                                  re =        mod_object$re, 
                                  vars =      mod_object$vars, 
                                  mod =       mod_object$mod , 
                                  alt_names = mod_object$alt_names, 
                                  theta =     mod_object$theta, 
                                  theta_0 =   mod_object$theta_0, 
                                  ll =        mod_object$ll, 
                                  H =         mod_object$H, 
                                  J =         mod_object$J, 
                                  grad =      mod_object$grad, 
                                  fit =       mod_object$fit, 
                                  vv =        mod_object$vv, 
                                  vv2 =       mod_object$vv2, 
                                  control =   mod_object$control, 
                                  info =      mod_object$info, 
                                  gradEL =    mod_object$gradEL, 
                                  HessEL =    mod_object$HessEL)
  
  if(!is.null(mod_object$mod)){
    out_mod <- mod_cl$new(Hb = as.matrix(mod_object$mod$Hb), 
                                    fb = as.matrix(mod_object$mod$fb), 
                                    HO = as.matrix(mod_object$mod$HO), 
                                    fO = as.matrix(mod_object$mod$fO), 
                                    HL = as.matrix(mod_object$mod$HL), 
                                    fL = as.matrix(mod_object$mod$fL), 
                                    alt = mod_object$mod$alt, 
                                    ordered = mod_object$mod$ordered)
    out$mod <- out_mod
  }
  
  if(!is.null(mod_object$control)){
    out_control <- control_cl$new(control_nlm = mod_object$control$control_nlm, 
                                            probit = mod_object$control$probit, 
                                            approx_method = mod_object$control$approx_method, 
                                            hess = as.logical(mod_object$control$hess), 
                                            pairs_list = mod_object$control$pairs_list, 
                                            el = as.logical(mod_object$control$el), 
                                            control_weights = mod_object$control$control_weights, 
                                            control_simulation = mod_object$control$control_simulation, 
                                            nCores = mod_object$control$nCores, 
                                            normalize = if(!is.null(mod_object$control$normalize)) as.logical(mod_object$control$normalize) else FALSE)
    
    out$control <- out_control
  }
  
  if(!is.null(mod_object$data_raw)){
    out_data_raw <- data_raw_cl$new(df = mod_object$data_raw, 
                                              alt_names = if(!is.null(mod_object$alt_names)) mod_object$alt_names else as.character(sort(unique(mod_object$data_raw$choice))), 
                                              id = if(!is.null(mod_object$info$setup_input$ids)) mod_object$info$setup_input$ids else "id_macml", 
                                              choice = if(!is.null(mod_object$info$setup_input$choice)) mod_object$info$setup_input$choice else "choice", 
                                              ordered = mod_object$mod$ordered)
    out$data_raw <- out_data_raw
  }
  
  return(out)
}
