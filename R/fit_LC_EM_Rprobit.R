#' Fit a Probit Model including latent classes. 
#' 
#' @description 
#' This function fits a probit model with a number of latent classes using the EM algorithm.
#' DOES NOT USE PARALLELIZATION!!!
#' 
#' @param Rprobit_obj 
#' An \code{\link{Rprobit_cl}} object.
#' @param init_method 
#' A \code{character}, setting the method for initialization of the likelihood 
#' optimization.
#' Can be one of 
#' * \code{"random"} (default) 
#' * or \code{"subset"}.
#' @param control_nlm 
#' Optionally a named \code{list} of control parameters passed on to 
#' \link[stats]{nlm}.
#' By default, \code{control_nlm = NULL}.
#' 
#' @param cml_pair_type 
#' An integer. 
#' cml_pair_type = 0 (full pairwise), = 1 (adjacent pairwise looped), = 2 (adjacent pairwise CML)
#' By default, \code{cml_pair_type = NULL}, such that the information in \code{Rprobit_obj} is used. 
#' @param control_EM
#' A list containing iter_lim (max number of iterations), iter_lim_one_step (max iterations in one step),
#' tol (tolerance for exiting)
#' 
#' @return 
#' A \code{\link{Rprobit_cl}} object.
#' 
#' @export

fit_LC_EM_Rprobit <- function(Rprobit_obj, init_method = "random", control_nlm = NULL, cml_pair_type = NULL,control_EM = NULL) {
  
  ## set controls for EM algorithm, if not supplied 
  if (is.null(control_EM)){
    control_EM$iter_lim <- 20 # total number of iterations
    control_EM$iter_lim_one_step <- 5 # number of iteration steps in one optimization
    control_EM$tol <- 0.00001 # loop ended, if criterion function changes less than tol. 
    control_EM$optim_all <- FALSE # (FALSE) (default: cycle over classes, optimising one per iteration)
                                  # (TRUE) optimize over all classes in each iteration
  }
  
  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  mod_orig <- Rprobit_o$mod$clone(deep=TRUE)
  data_raw_orig = Rprobit_o$data_raw$clone(deep=TRUE)
  ### check if there are random effects. If not, no need to go to complicated calculations, remove panel structure.
  if ((Rprobit_o$mod$ordered == FALSE)&&(is.null(Rprobit_o$re)==TRUE)){
    if (!is.null(Rprobit_o$data_raw)){
      df = Rprobit_o$data_raw$df
      df[,Rprobit_o$data_raw$id] = c(1:dim(Rprobit_o$data_raw$df)[1])
      Rprobit_o$data_raw$set_df(df)
      N <- dim(Rprobit_o$data_raw$df)[1]
      # Rprobit_o$mod$N  = N
      # Rprobit_o$mod$Tp = sample(1,N,replace=TRUE)
    }
    if (!is.null(Rprobit_o$data)){
      data = list()
      cur = 0
      for (ja in 1:length(Rprobit_o$data$data)){
        data_n = Rprobit_o$data$data[[ja]]
        for (jb in 1:length(data_n$y)){
          dat = list()
          dat$X[[1]] = data_n$X[[jb]]
          dat$y = data_n$y[jb]
          cur = cur+1
          data[[cur]] <- dat
        }
      }
      Rprobit_o$data$set_data(data)
    }

  }

  #####################################################
  ## set up helper functions for optimization.    #####
  #####################################################
  # helper functions needed in the algorithm 
  # calculate conditional probabilities.
  col_cond_prob <- function(perc,pi){
    N = dim(perc)[1]
    K = dim(perc)[2]
    perc <- perc * (matrix(1,N,1) %*% pi)
    sumperc= apply(perc,1,sum)
    for (k in 1:K){
      perc[,k]=perc[,k]/sumperc
    }
    return( perc)
  }
  
  # evaluate negative log-likelihood of latent class model  
  eval_LC_ll <- function(ll_s,pi){
    pr_s <- exp(-ll_s)
    LN <- sum( log(pr_s %*% pi))
    return ((-1)*LN)
  }
  
  #######################################################
  ### The main EM function                          #####
  #######################################################
  EM_loop <- function(theta,K,param_one,data_tr,mod,control,control_EM){
    
    N = data_tr$N
    # start EM computations 
    # steps:
    
    # calculate conditional class membership probabilities 
    ll_s <- matrix(0,N,K)
    for (kchoose in 1:K){
      outk <- ll_macml_marginal(theta = theta[(kchoose-1)*param_one+c(1:param_one)], data_obj =data_tr, mod = mod, control = control)
      ll_s[,kchoose] <- unlist(attr(outk,"llEL"))
    }
    
    # calculate conditional distribution
    pi_param <- c(0,theta[K*param_one +c(1:(K-1))])
    pi_est<- exp(pi_param)/sum(exp(pi_param))
    perc_est <- col_cond_prob(exp(-ll_s),pi_est)
    pi_est <- apply(perc_est,2,sum)/N
    
    # initial value of negative log-likelihood 
    ll_old <- eval_LC_ll(ll_s,pi_est)
    print.level <- control$control_nlm$print.level 
    
     if (print.level>0){
      cat('iteration = 0',"\n")
      cat('Parameter: ',"\n")
      cat(theta,"\n")
      cat('Function Value',"\n")
      cat(ll_old,"\n")
    }
   
    tol <- control_EM$tol

    ll <- ll_old 
    ll_old <- ll_old+ 1 
    it <- 0 
    iter_lim <- control_EM$iter_lim
    iter_lim_one_step <- control_EM$iter_lim_one_step
    
    while ((ll<ll_old-tol)&&(it<iter_lim)){
      it <- it+1 
      
      theta_old <- theta
      # check convergence criterion 
      ll_old <- ll 
      # start EM computations 
      perc_est <- col_cond_prob(exp(-ll_s),pi_est)
      
      pi_est <- apply(perc_est,2,sum)/N
      # maximize the full data likelihood 
      eval_ll <- function(theta,k,ll_s,perc,data_tr, mod, control){
        out <- ll_macml_marginal(theta, data_obj =data_tr, mod = mod, control = control)
        ll_s[,k] <- attr(out,"llEL")
        
        ll_new = sum(ll_s*perc)
        # calculate gradient 
        gr <- matrix(unlist(attr(out,"gradEL")),nrow= length(theta))
        attr(ll_new,"gradient") <- gr %*% perc[,k]
        attr(ll_new,"ll_s") <- ll_s
        return(ll_new)
      }
      
      # two different approaches: (A) optimize for all classes (B) Optimize only one class per iteration.
      if (control_EM$optim_all == TRUE){
        
        for (kchoose in 1:K){
            outnlm <- stats::nlm(f  = eval_ll,
                             p = theta[(kchoose-1)*param_one+c(1:param_one)],
                             k = kchoose, 
                             ll_s = ll_s,
                             perc = perc_est, 
                             data_tr   = data_tr,
                             mod       = mod,
                             control   = control,
                             print.level = 0,
                             iterlim     = iter_lim_one_step)
        
            # update estimate for ll_s
            out <- ll_macml_marginal(outnlm$estimate, data_obj =data_tr, mod = mod, control = control)
            ll_s[,kchoose] <- attr(out,"llEL")
            theta[(kchoose-1)*param_one+c(1:param_one)] <- outnlm$estimate
        } 
      } else {
        
          # draw class to maximize 
          kchoose <- kchoose+1
          if (kchoose == K+1){kchoose <- 1}
          
          outnlm <- stats::nlm(f  = eval_ll,
                        p = theta[(kchoose-1)*param_one+c(1:param_one)],
                        k = kchoose, 
                        ll_s = ll_s,
                        perc = perc_est, 
                        data_tr   = data_tr,
                        mod       = mod,
                        control   = control,
                        print.level = 0,
                        iterlim     = iter_lim_one_step)
          
          # update estimate for ll_s
          out <- ll_macml_marginal(outnlm$estimate, data_obj =data_tr, mod = mod, control = control)
          ll_s[,kchoose] <- attr(out,"llEL")
          theta[(kchoose-1)*param_one+c(1:param_one)] <- outnlm$estimate
      }
      
      
      # calculate conditional class membership probabilities 
      perc_est <- col_cond_prob(exp(-ll_s),pi_est)
      
      # update estimate membership probabilities
      pi_est <- apply(perc_est,2,sum)/N
      pi_est_norm <- pi_est/pi_est[1]
      gamma_est <- log(pi_est_norm[-1])
      
      theta[K*param_one + c(1:length(gamma_est))] <- gamma_est
      
      # evaluate new value of criterion function 
      ll <- eval_LC_ll(ll_s,pi_est)
      if (print.level>0){
        cat(sprintf('iteration = %d',it),"\n")
        cat('Parameter: ',"\n")
        cat(theta,"\n")
        cat('Function Value',"\n")
        cat(ll, "\n")
        cat('Update Step',"\n")
        cat(theta-theta_old,"\n\n")
      }
    } 
    
    # prepare output of function 
    # convert pi_est in parameters. 
    pi_est <- pi_est/pi_est[1]
    gamma_est <- log(pi_est[-1])
    
    theta[K*param_one + c(1:length(gamma_est))] <- gamma_est
    
    return (list(theta_est =theta, it = it))
  }
  
  

  ### extract print level
  print.level           <- 0

  if(!is.null(control_nlm$print.level)){
    print.level    <- control_nlm$print.level
  } else if(!is.null(Rprobit_o$control$control_nlm$print.level)){
    print.level    <- Rprobit_o$control$control_nlm$print.level
  }

  ### make sure that 'Rprobit_o' is of class 'Rprobit'
  if(!inherits(Rprobit_o, "Rprobit_cl")) {
    stop("'Rprobit_o' is not a Rprobit_cl object")
  }

  ### unpack the control variables
  probit            <- FALSE
  approx_method     <- "SJ"
  hessian           <- FALSE

  normalize         <- FALSE # if not supplied, default is not to normalize.
  if(!is.null(Rprobit_o$control$approx_method))       approx_method     <- Rprobit_o$control$approx_method
  if(!is.null(Rprobit_o$control$probit))              probit            <- Rprobit_o$control$probit
  if(!is.null(Rprobit_o$control$hess))                hessian           <- Rprobit_o$control$hess
  if(!is.null(Rprobit_o$control$normalize))           normalize         <- Rprobit_o$control$normalize
  Rprobit_o$control$approx_method <- approx_method
  Rprobit_o$control$probit        <- probit
  Rprobit_o$control$normalize     <- normalize

  ### Check if parallel framework should be used and transform data in that case
  use_parallel <- FALSE
  if (!is.null(Rprobit_o$control$nCores)) {
    if (Rprobit_o$control$nCores > 1) {
      use_parallel <- TRUE
    }
  }


  ### select likelihood function
  if(probit)  ll_function <- ll_probit

  if(!probit){
    if (normalize == TRUE){
      ll_function <- ll_macml_norm
    } else {
      ll_function <- ll_macml_marginal
    }
  }

  if(Rprobit_o$mod$ordered == TRUE) ll_function <- ll_macml_o


  ### check if list version of data is available:
  if(is.null(Rprobit_o[["data"]]) & !is.null(Rprobit_o[["data_raw"]])){
    data_from_data_raw_model <- function(Rprobit_o){
      read_formula_out  <- read_formula(Rprobit_o$form)
      ASC               <- read_formula_out$ASC
      allvars           <- read_formula_out$allvars
      choice            <- all.vars(Rprobit_o$form)[1]
      norm_alt          <- Rprobit_o$info$setup_input$norm_alt
      alt               <- Rprobit_o$mod$alt

      data_raw_to_data_out = data_raw_to_data(data_raw = Rprobit_o$data_raw, allvars = allvars, choice = choice, re = Rprobit_o$re, norm_alt = norm_alt, alt = alt, ASC = ASC)
      #data_raw_to_data_out  <- data_raw_to_data(data_raw = Rprobit_o$data_raw, allvars = allvars, choice = choice, re = Rprobit_o$re, norm_alt = norm_alt, alt = alt, ASC = ASC)

      return(data_raw_to_data_out)
    }
    if (Rprobit_o$mod$ordered){
      data_raw_obj = Rprobit_o$data_raw$clone()
      data = data_cl$new(ordered = TRUE,vars = data_raw_obj$dec_char)
      ids = data_raw_obj$df[,data_raw_obj$id]
      unique_ids = unique(ids)
      data_df = list()
      for (j in 1:length(unique_ids)){
        ind = which(ids == unique_ids[j])
        data_df[[j]] = list(X = data.matrix(data_raw_obj$df[ind,data_raw_obj$dec_char]), y = data.matrix(data_raw_obj$df[ind,data_raw_obj$choice]))
      }
      data$set_data(data_df)
      Rprobit_o$data <- data
    } else {
      Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
    }
  }


  ### transform data for efficient likelihood computation
  if (Rprobit_o$mod$ordered == FALSE) {
    data_tr <- Rprobit_o$data$clone()
    data_tr$data <- substract_choice_regressor_from_data(Rprobit_o$data$data)
  } else {
    data_tr <- Rprobit_o$data # for ordered probit there is no need to difference with respect to chosen alternative.
  }


  ### extract nlm controls
  {
    typsize_nlm <- rep(1, length(Rprobit_o$theta))

    fscale_nlm            <- 1
    if(!is.null(control_nlm$fscale)){
      fscale_nlm <- control_nlm$fscale
    } else if (!is.null(Rprobit_o$control$control_nlm$fscale)) {
      fscale_nlm <- Rprobit_o$control$control_nlm$fscale
    }

    print.level_nlm       <- min(print.level, 2)

    ndigit_nlm            <- 12
    if(!is.null(control_nlm$ndigit)){
      ndigit_nlm <- control_nlm$ndigit
    } else if (!is.null(Rprobit_o$control$control_nlm$ndigit)) {
      ndigit_nlm <- Rprobit_o$control$control_nlm$ndigit
    }
    gradtol_nlm           <- 1e-06
    if(!is.null(control_nlm$gradtol)){
      gradtol_nlm <- control_nlm$gradtol
    } else if (!is.null(Rprobit_o$control$control_nlm$gradtol)) {
      gradtol_nlm <- Rprobit_o$control$control_nlm$gradtol
    }

    stepmax_nlm <- max(1000 * sqrt(sum((Rprobit_o$theta / typsize_nlm)^2)), 1000)
    if (is.na(stepmax_nlm)) {
      stepmax_nlm <- 1000
    }
    if (!is.null(control_nlm$stepmax)) {
      stepmax_nlm <- control_nlm$stepmax
    } else if (!is.null(Rprobit_o$control$control_nlm$stepmax)) {
      stepmax_nlm <- Rprobit_o$control$control_nlm$stepmax
    }

    steptol_nlm <- 1e-06
    if (!is.null(control_nlm$steptol)) {
      steptol_nlm <- control_nlm$steptol
    } else if (!is.null(Rprobit_o$control$control_nlm$steptol)) {
      steptol_nlm <- Rprobit_o$control$control_nlm$steptol
    }

    iterlim_nlm <- 1000
    if (!is.null(control_nlm$iterlim)) {
      iterlim_nlm <- control_nlm$iterlim
    } else if (!is.null(Rprobit_o$control$control_nlm$iterlim)) {
      iterlim_nlm <- Rprobit_o$control$control_nlm$iterlim
    }

    check.analyticals_nlm <- FALSE
    if (!is.null(control_nlm$check.analyticals)) {
      check.analyticals_nlm <- control_nlm$check.analyticals
    } else if (!is.null(Rprobit_o$control$control_nlm$check.analyticals)) {
      check.analyticals_nlm <- Rprobit_o$control$control_nlm$check.analyticals
    }

    force_hess <- FALSE
    if (!is.null(control_nlm$force_hess)) {
      force_hess <- control_nlm$force_hess
    } else if (!is.null(Rprobit_o$control$control_nlm$force_hess)) {
      force_hess <- Rprobit_o$control$control_nlm$force_hess
    }


    Rprobit_o$control$control_nlm <- list(
      typsize = typsize_nlm,
      fscale = fscale_nlm,
      print.level = print.level_nlm,
      ndigit = ndigit_nlm,
      gradtol = gradtol_nlm,
      stepmax = stepmax_nlm,
      steptol = steptol_nlm,
      iterlim = iterlim_nlm,
      check.analyticals = check.analyticals_nlm,
      force_hess = force_hess
    )
  }


  ### only optimise when iterlim > 0
  if (Rprobit_o$control$control_nlm$iterlim > 0) {

    ### create initial parameter vector
    if (print.level > 0) {
      if (!is.null(Rprobit_o$theta)) {
        print("Check initial parameter values")
      } else {
        print("Create initial parameter values")
      }
    }


    ### TODO: update initialization 
    ### only random initialization for now !! Rprobit_o$theta <- init_Rprobit(ll_function = ll_function, Rprobit_obj = Rprobit_o, init_method = init_method, data_tr = data_tr)
    ### check whether the dimension of theta matches the dimensions given in mod
    lth = Rprobit_o$mod$lthb + Rprobit_o$mod$lthO + Rprobit_o$mod$lthL
    if (!is.null(Rprobit_o$mod$num_class)){
      K <- Rprobit_o$mod$num_class
    } else {
      cat('No number of classes given, assuming K=2!')
      K <- 2
    }
    
    tot_params <- (lth+1)*K -1
    
    Rprobit_o$theta <- rnorm(tot_params)

    if (Rprobit_o$mod$ordered == TRUE){
      lth = lth + Rprobit_o$mod$alt-1 # add parameters for tauk.
    }
    
    Rprobit_o$theta <- rnorm(tot_params)
    
    if (init_method == "kmeans") {
      Rprobit_o$theta <- init_LC_Rprobit(Rprobit_o,mod=Rprobit_o$mod,control= Rprobit_o$control)
    }
    
    if (length(Rprobit_o$theta) != tot_params){
      stop(paste0("Number of parameters does not match the structure of the mod object. "))
    }

    ### perform numerical optimization and save results
    Rprobit_o$control$hess <- FALSE
    if (print.level > 0) {
      print("Start optimisation...")
    }
    time_1 <- Sys.time()
    ### calculate weights, if necessary
    
    if (!is.null(cml_pair_type)){
      # add CML type
      Rprobit_o$control$control_weights$cml_pair_type = cml_pair_type
    }
    if (!is.null(Rprobit_o$control$control_weights)){
      Rprobit_o       <- CML_weights(Rprobit_o, control_weights = Rprobit_o$control$control_weights)
    }
    
    Rprobit_o$data = data_tr
    
    N <- data_tr$N
    
    
    out_opt <- EM_loop(Rprobit_o$theta,K,param_one=lth,data_tr,Rprobit_o$mod,Rprobit_o$control,control_EM)
    time_2          <- Sys.time()
    
    # calculate result at estimate
    Rprobit_o$control$hess <- FALSE
    Rprobit_o$theta <- out_opt$theta_est
    neg_ll_fit <- ll_macml_LC(Rprobit_o$theta,data_obj = Rprobit_o$data, mod=Rprobit_o$mod,control = Rprobit_o$control)
    Rprobit_o$ll <- (-1) * neg_ll_fit
    #Rprobit_o$grad <- attr(neg_ll_fit, "gradient")
    Rprobit_o$control$hess <- hessian

    Rprobit_o$fit <- approx_method
    Rprobit_o$control$hess <- hessian
    if (!is.null(Rprobit_o$info)) {
      Rprobit_o$info$estimation_time <- difftime(time_2, time_1, units = "secs")
      Rprobit_o$info$nlm_info <- out_opt$theta_est
    } else {
      Rprobit_o$info <- list(
        estimation_time = difftime(time_2, time_1, units = "secs"),
        nlm_info <- out_opt$theta_est
      )
    }

    if (out_opt$it == control_EM$iter_lim){ # iterated till the bitter end
      warning("Used up all iterations, solution probably not found.")
      }
    ### print nlm exit codes
    #if (nlm_out$code == 2) {
    #  warning("Exit code 2: successive iterates within tolerance, current iterate is probably solution")
    #} else if (nlm_out$code == 3) {
    #  warning("Exit code 3: last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.")
    #} else if (nlm_out$code == 4) {
    #  warning("Exit code 4: iteration limit exceeded.")
    #} else if (nlm_out$code == 5) {
    #  warning("Exit code 5: maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.")
    #}

    ### don't compute Hessian Matrix if nlm_code suggests estimation was not successful
    #if (nlm_out$code > 3) {
    #  if (!(Rprobit_o$control$control_nlm$force_hess)) {
    #    hessian <- FALSE
    #  }
    #}
  } else {
    Rprobit_o$control$hess <- FALSE
    neg_ll_fit <- ll_macml_LC(Rprobit_o$theta,data_obj = Rprobit_o$data, mod=Rprobit_o$mod,control = Rprobit_o$control)
    Rprobit_o$ll <- (-1) * neg_ll_fit
    #Rprobit_o$grad <- attr(neg_ll_fit, "gradient")
    Rprobit_o$control$hess <- hessian
  }


  ### computation of variance

  if (print.level > 0) {
    print("Estimate Hessian matrix...")
  }

  ### register parallel setup if necessary
  #re_register_dopar(Rprobit_o = Rprobit_o, use_parallel = use_parallel)

  time_1                    <- Sys.time()
  #if(Rprobit_o$control$approx_method!="TVBS"){

  if (Rprobit_o$control$hess == TRUE){
    Rprobit_o$H               <- attr(neg_ll_fit, "hessian1")
  } else {
    Rprobit_o$H <- attr(neg_ll_fit, "hessian1")
  }
  if (Rprobit_o$control$probit == TRUE){
    Rprobit_o$J               <- attr(neg_ll_fit, "hessian2")
  } else {
    Rprobit_o$J               <- attr(neg_ll_fit, "hessian2")
  }

  ### save individual contributions to gradient and hessian
  #if ("gradEL" %in% names(attributes(neg_ll_fit))) {
  #  Rprobit_o$gradEL <- attr(neg_ll_fit, "gradEL")
  #}
  #if ("HessEL" %in% names(attributes(neg_ll_fit))) {
  #  Rprobit_o$HessEL <- attr(neg_ll_fit, "HessEL")
  #}

  ### if covariance matrix H is not pos.definite, adjust
  H_fixed                   <- adjust_variance_matrix(Rprobit_o$H)
  if (H_fixed$ok == FALSE){
    message("Calculated variance matrix H needed to be adjusted due to negative eigenvalue(s)!")
  }
  inverted_H                <- H_fixed$mat_inv

  ### if covariance matrix J is not pos.definite, adjust
  S = adjust_variance_matrix(Rprobit_o$J)
  if (S$ok == FALSE){
    message("Calculated variance matrix J needed to be adjusted due to eigenvalue(s) close to zero!")
  }
  if (Rprobit_o$control$probit == TRUE) {
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

  if (print.level > 0) {
    print("...done.")
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
