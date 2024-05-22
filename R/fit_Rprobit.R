#' Fit a Probit Model
#' 
#' @description 
#' This function fits a probit model.
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
#' @return 
#' A \code{\link{Rprobit_cl}} object.
#' 
#' @export

fit_Rprobit <- function(Rprobit_obj, init_method = "random", control_nlm = NULL, cml_pair_type = NULL) {
  
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
      # Rprobit_o$mod$N = length(data)
      # Rprobit_o$mod$Tp = sample(1,length(data),replace=TRUE)
    }

  }

  ### if ordered, then probit estimation does not work, we always use CML!
  if (Rprobit_o$mod$ordered == TRUE){
    Rprobit_o$control$probit = FALSE
  }
  ### check model first ####
  # ok = check_mod_par(Rprobit_o$mod)

  ### extract print level
  print.level           <- 0

  if(!is.null(control_nlm$print.level)){
    print.level    <- control_nlm$print.level
  } else if(!is.null(Rprobit_o$control$control_nlm$print.level)){
    print.level    <- Rprobit_o$control$control_nlm$print.level
  }



  ### define optimizer function
  optimizer <- function(ll_function, Rprobit_o, data){


    ### define nlm parameters
    p_nlm <- Rprobit_o$theta


    ### initiate sinking the nlm output if print.level>2
    if (print.level > 2) {
      sink_name <- paste(paste(unlist(strsplit(x = Rprobit_o$info$name, split = " ")), collapse = "_"), "_nlm_output_", Sys.Date(), ".txt", sep = "")
      sink(
        file = sink_name,
        split = TRUE,
        append = TRUE
      )
    }
    Rprobit_o$control$control_nlm$typsize = rep(1,length(p_nlm))

    ### optimise
    out <- stats::nlm(f                  = ll_function,
                      p                  = p_nlm,
                      data               = data,
                      mod                = Rprobit_o$mod,
                      control            = Rprobit_o$control,
                      typsize            = Rprobit_o$control$control_nlm$typsize,
                      fscale             = Rprobit_o$control$control_nlm$fscale,
                      print.level        = min(Rprobit_o$control$control_nlm$print.level, 2),
                      ndigit             = Rprobit_o$control$control_nlm$ndigit,
                      gradtol            = Rprobit_o$control$control_nlm$gradtol,
                      stepmax            = Rprobit_o$control$control_nlm$stepmax,
                      steptol            = Rprobit_o$control$control_nlm$steptol,
                      iterlim            = Rprobit_o$control$control_nlm$iterlim,
                      check.analyticals  = Rprobit_o$control$control_nlm$check.analyticals)

    if(print.level>2){
      ### if print.level > 2, stop saving the text output
      sink(file = NULL)
      if (print.level == 3) {
        ### if print.level==3, convert text output to readable data.frame and save it in a .csv file
        parameter_names <- NULL
        if (Rprobit_o$mod$lthb > 0) {
          if (length(Rprobit_o$vars) == Rprobit_o$mod$lthb) {
            parameter_names <- c(parameter_names, Rprobit_o$vars)
          } else {
            parameter_names <- c(parameter_names, paste("b_", 1:Rprobit_o$mod$lthb, sep = ""))
          }
        }
        if(Rprobit_o$mod$lthO>0) parameter_names <- c(parameter_names, paste("omega_", formatC(1:Rprobit_o$mod$lthO, width=ceiling(log10(Rprobit_o$mod$lthO)), flag="0"), sep = ""))
        if(Rprobit_o$mod$lthL>0) parameter_names <- c(parameter_names, paste("sigma_", formatC(1:Rprobit_o$mod$lthL, width=ceiling(log10(Rprobit_o$mod$lthL)), flag="0"), sep = ""))


        nlm_out_2_data_frame(file_name = sink_name,
                             save_in_file = paste(strsplit(x = sink_name, split = "_nlm_output_")[[1]][1],
                                                  strsplit(x = strsplit(x = sink_name, split = "_nlm_output_")[[1]][2], split = ".txt")[[1]],
                                                  sep = "_Parameters_"
                             ),
                             return_dataframe = FALSE,
                             parameter_names = parameter_names
        )
      }
      ### remove the .txt file
      unlink(x = sink_name)
    }

    return(out)
  } # end optimizer function

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
      ll_function <- ll_macml
    }
  }

  # for ordered probit, use different CML implementations 
  if(Rprobit_o$mod$ordered == TRUE){ # if error includes state space autocorrelation 
    if (inherits(mod_orig,"mod_StSp_cl")){ 
      ll_function <- ll_macml_o_StSp
    } else {  # else no autocorrelation. 
      if (inherits(mod_orig,"mod_AR_cl")){ 
        ll_function <- ll_macml_o_AR
      } else {
        ll_function <- ll_macml_o
      }
    }
  } 


  ### define functions for parallel setup
  ### to clean up parallel computing:
  #unregister_dopar <- function() {
  #  env <- foreach:::.foreachGlobals
  #  rm(list=ls(name=env), pos=env)
  # }
  #
  # re_register_dopar <- function(Rprobit_o, use_parallel){
  #  if(use_parallel){
  #    unregister_dopar()
  #
  #    clusters <- parallel::makeCluster(Rprobit_o$control$nCores)
  #    doSNOW::registerDoSNOW(clusters)
  #  }
  # }
  #
  # stop_parallel <- function(Rprobit_o, use_parallel){
  #  if(use_parallel){
  #    unregister_dopar()
  #    clusters <- parallel::makeCluster(Rprobit_o$control$nCores)
  #    parallel::stopCluster(clusters)
  #  }
  #}

  ### build parallel setup
  if(use_parallel){

    ### check number of available cores
    Rprobit_o$control$nCores <- min(parallel::detectCores()-1, Rprobit_o$control$nCores)

    ### need dopar function from foreach package
    `%dopar%` <- foreach::`%dopar%`


    ### function to create batches
    Rprobit_data_batches <- function(Tp, n_Batch = 1) {
      out <- list()

      if(n_Batch>1){

        ### order of Tps
        order_Tp <- order(Tp, decreasing = TRUE)

        ### total number of observations per batch
        n_Tp <- numeric(n_Batch)

        ### initialize list
        for (i in 1:n_Batch) {
          out[[i]] <- numeric(0)
        }

        ### sort each person to one batch
        for (i in 1:length(Tp)) {
          add_to <- which.min(n_Tp)

          out[[add_to]] <- c(out[[add_to]], order_Tp[i])
          n_Tp[add_to] <- n_Tp[add_to] + Tp[order_Tp[i]]
        }
      } else {
        out[[1]] <- 1:length(Tp)
      }
      return(out)
    }

    ### create batches
    #macml_batches <- Rprobit_data_batches(Tp = Rprobit_o$mod$Tp, n_Batch = Rprobit_o$control$nCores)

    ll_Rpar <- function(theta, data, mod, control){
      ### initialise output
      ll_par <- 0

      ### create batches
      Rprobit_batches_parallel_int <- Rprobit_data_batches(Tp = data$Tp, n_Batch = control$nCores)

      ### internal functions work with only one core
      save_cores      <- control$nCores
      control$nCores  <- 1

      ### parallelise logLik code
      clusters <- parallel::makeCluster(save_cores)
      doSNOW::registerDoSNOW(clusters)

      data_tr_batch = data$clone() # clone data to specify different data batches
      ll_list_par <- foreach::foreach(i_batch=1:length(Rprobit_batches_parallel_int)) %dopar% {
        data_tr_batch$set_data(data$data[Rprobit_batches_parallel_int[[i_batch]]])
        if (R6::is.R6(mod)) {
          mod_batch <- mod$clone() ### TODO: does not work with manually defined mod objects!
        } else {
          mod_batch <- mod
        }
        # mod_batch$Tp  <- mod_batch$Tp[Rprobit_batches_parallel_int[[i_batch]]]
        # mod_batch$N   <- length(Rprobit_batches_parallel_int[[i_batch]])
        control_batch <- control$clone()
        if (length(control$pairs_list) > 0) {
          control_batch$pairs_list <- control$pairs_list[Rprobit_batches_parallel_int[[i_batch]]]
        }
        if (control$probit) {
          ll_rpar_out <- Rprobit::ll_probit(theta = theta, data_obj = data_tr_batch, mod = mod_batch, control = control_batch)
        } else {
          ll_rpar_out <- Rprobit::ll_macml(theta = theta, data_obj = data_tr_batch, mod = mod_batch, control = control_batch)
        }
        ll_rpar_out
      }

      ### stop clusters
      parallel::stopCluster(clusters)

      ### initialise attributes
      attr_ll <- names(attributes(ll_list_par[[1]]))
      for (i_attr in 1:length(attr_ll)) {
        if (!(attr_ll[i_attr] %in% c("gradEL", "HessEL"))) {
          attr(ll_par, attr_ll[i_attr]) <- attr(ll_list_par[[1]], attr_ll[i_attr]) * 0
        }
      }

      el <- FALSE
      if (!is.null(Rprobit_o$control$el)) {
        el <- Rprobit_o$control$el
      }
      if (el == TRUE | el == 1) {
        attr(ll_par, "gradEL") <- list()
        #attr(ll_par, "HessEL") <- list()
      }

      ###
      for(i_batch in 1:length(Rprobit_batches_parallel_int)){
        ll_par <- ll_par + ll_list_par[[i_batch]]

        for(i_attr in 1:length(attr_ll)){
          if(!(attr_ll[i_attr] %in% c("gradEL", "HessEL"))){
            attr(ll_par, attr_ll[i_attr]) <- attr(ll_par, attr_ll[i_attr]) + attr(ll_list_par[[i_batch]], attr_ll[i_attr])
          }
        }
        if ((el == TRUE | el == 1) & "gradEL" %in% attr_ll) {
          for (i_attr_n in 1:length(attr(ll_list_par[[i_batch]], "gradEL"))) {
            attr(ll_par, "gradEL")[[Rprobit_batches_parallel_int[[i_batch]][i_attr_n]]] <- attr(ll_list_par[[i_batch]], "gradEL")[[i_attr_n]]
          }
        }
        if ((el == TRUE | el == 1) & "HessEL" %in% attr_ll) {
          for (i_attr_n in 1:length(attr(ll_list_par[[i_batch]], "gradEL"))) {
            attr(ll_par, "HessEL")[[Rprobit_batches_parallel_int[[i_batch]][i_attr_n]]] <- attr(ll_list_par[[i_batch]], "HessEL")[[i_attr_n]]
          }
        }
      }
      control$nCores <- save_cores

      return(ll_par)
    }


    ll_function <- ll_Rpar
  } # end if use parallel


  ### if no mixed effects, set HO and fO to zero matrices
  #if (Rprobit_o$mod$lthO == 0) {
  #  Rprobit_o$mod$HO <- matrix(0, 0, 0)
  #  Rprobit_o$mod$fO <- matrix(0, 0, 0)
  #}

  ### if no mixed effects, set HO and fO to zero matrices
  #if (any(is.na(Rprobit_o$mod$HO))) {
  #  Rprobit_o$mod$HO <- matrix(0, 0, 0)
  #}
  #if (any(is.na(Rprobit_o$mod$fO))) {
  #  Rprobit_o$mod$fO <- matrix(0, 0, 0)
  #}


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
    # no distinction necessary any more. 
    #if (Rprobit_o$mod$ordered){
    #  data_raw_obj = Rprobit_o$data_raw$clone()
    #  data = data_cl$new(ordered = TRUE,vars = data_raw_obj$dec_char)
    #  ids = data_raw_obj$df[,data_raw_obj$id]
    #  unique_ids = unique(ids)
    #  data_df = list()
    #  for (j in 1:length(unique_ids)){
    #    ind = which(ids == unique_ids[j])
    #    data_df[[j]] = list(X = data.matrix(data_raw_obj$df[ind,data_raw_obj$dec_char]), y = data.matrix(data_raw_obj$df[ind,data_raw_obj$choice]))
    #  }
    #  #vars = allvars[[1]]
    #  #alt_names = NULL
    #  data$set_data(data_df)
    #  Rprobit_o$data <- data
    #} else {
      Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
    #}
  }


  ### transform data for efficient likelihood computation
  if (Rprobit_o$mod$ordered == FALSE) {
    data_tr <- Rprobit_o$data$clone()
    data_tr$data <- substract_choice_regressor_from_data(Rprobit_o$data$data)
  } else {
    data_tr <- Rprobit_o$data$clone() # for ordered probit there is no need to difference with respect to chosen alternative.
  }


  ### extract nlm controls
  {
    typsize_nlm <- rep(1, Rprobit_o$mod$lthb + Rprobit_o$mod$lthO + Rprobit_o$mod$lthL) # rep(1, length(Rprobit_o$theta))
    #if (!is.null(control_nlm$typsize)) {
    #  typsize_nlm <- control_nlm$typsize
    #} else if (!is.null(Rprobit_o$control$control_nlm$typsize)) {
    #  typsize_nlm <- Rprobit_o$control$control_nlm$typsize
    #}

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

    #re_register_dopar(Rprobit_o = Rprobit_o, use_parallel = use_parallel)
    # always initialize randomly
    lth <- Rprobit_o$mod$lthb + Rprobit_o$mod$lthO + Rprobit_o$mod$lthL
    #Rprobit_o$theta <- rnorm(lth)
    Rprobit_o$theta <- init_Rprobit(ll_function = ll_function, Rprobit_obj = Rprobit_o, init_method = init_method, data_tr = data_tr)
    ### check whether the dimension of theta matches the dimensions given in mod
    lth = Rprobit_o$mod$lthb + Rprobit_o$mod$lthO + Rprobit_o$mod$lthL
    if (Rprobit_o$mod$ordered == TRUE){
      lth = lth + Rprobit_o$mod$alt-1 # add parameters for tauk.
      if (inherits(Rprobit_o$mod,"mod_StSp_cl")){
        s = -0.5 + sqrt(2*dim(Rprobit_o$mod$fL)[1]+.25)
        lth = lth+ 2*Rprobit_o$mod$dim_state*s
      }
      if (inherits(Rprobit_o$mod,"mod_AR_cl")){
        lth = lth+ Rprobit_o$mod$lag_length
      }
    }
    if (length(Rprobit_o$theta) != lth){
      stop(paste0("Number of parameters does not match the structure of the mod object."))
    }

    ### perform numerical optimization and save results
    Rprobit_o$control$hess <- FALSE
    if (print.level > 0) {
      print("Start optimisation...")
    }
    time_1 <- Sys.time()
    ### calculate weights, if necessary
    #re_register_dopar(Rprobit_o = Rprobit_o, use_parallel = use_parallel)
    
    if (!is.null(cml_pair_type)){
      # add CML type
      Rprobit_o$control$control_weights$cml_pair_type = cml_pair_type
    }
    if (!is.null(Rprobit_o$control$control_weights)){
      Rprobit_o       <- CML_weights(Rprobit_o, control_weights = Rprobit_o$control$control_weights)
    }
    
    ### register parallel setup if necessary
    #re_register_dopar(Rprobit_o = Rprobit_o, use_parallel = use_parallel)

    Rprobit_o$data = data_tr
    
    nlm_out         <- optimizer(ll_function = ll_function, Rprobit_o = Rprobit_o, data = Rprobit_o$data)
    time_2          <- Sys.time()
    Rprobit_o$ll    <- (-1)*(nlm_out$minimum)
    Rprobit_o$theta <- nlm_out$estimate
    Rprobit_o$grad <- nlm_out$gradient
    Rprobit_o$fit <- approx_method
    Rprobit_o$control$hess <- hessian
    if (!is.null(Rprobit_o$info)) {
      Rprobit_o$info$estimation_time <- difftime(time_2, time_1, units = "secs")
      Rprobit_o$info$nlm_info <- nlm_out
    } else {
      Rprobit_o$info <- list(
        estimation_time = difftime(time_2, time_1, units = "secs"),
        nlm_info = nlm_out
      )
    }

    ### print nlm exit codes
    if (nlm_out$code == 2) {
      warning("Exit code 2: successive iterates within tolerance, current iterate is probably solution")
    } else if (nlm_out$code == 3) {
      warning("Exit code 3: last global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function or steptol is too small.")
    } else if (nlm_out$code == 4) {
      warning("Exit code 4: iteration limit exceeded.")
    } else if (nlm_out$code == 5) {
      warning("Exit code 5: maximum step size stepmax exceeded five consecutive times. Either the function is unbounded below, becomes asymptotic to a finite value from above in some direction or stepmax is too small.")
    }

    ### don't compute Hessian Matrix if nlm_code suggests estimation was not successful
    if (nlm_out$code > 3) {
      if (!(Rprobit_o$control$control_nlm$force_hess)) {
        hessian <- FALSE
      }
    }
  } else {
    Rprobit_o$control$hess <- FALSE
    neg_ll_fit <- ll_function(theta = Rprobit_o$theta, data = data_tr, mod = Rprobit_o$mod, control = Rprobit_o$control)
    Rprobit_o$ll <- (-1) * neg_ll_fit
    Rprobit_o$grad <- attr(neg_ll_fit, "gradient")
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

  neg_ll_fit                <- ll_function(theta = Rprobit_o$theta, data = data_tr, mod = Rprobit_o$mod, control = Rprobit_o$control)
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

    #} else {
    #  Rprobit_o$control$hess    <- 1
    #  neg_ll_fit                <- ll_function(theta = Rprobit_o$theta, data = data_tr, mod = Rprobit_o$mod, control = Rprobit_o$control)
    #  Rprobit_o$J               <- attr(neg_ll_fit,"hessian1")
    #
    #  ### if TVBS is used, exact Hessian has to be computed numerically from the gradient
    #  ll_grad <- function(theta, data, mod, control){
    #    control$hess              <- 0
    #    ll_out                    <- ll_macml(theta = theta, data = data, mod = mod, control = control)
    #    out                       <- attr(ll_out, "gradient")
    #    return(out)
    #  }
    #  num_Hessian               <- numDeriv::jacobian(ll_grad, x = Rprobit_o$theta, data = data_tr, mod = Rprobit_o$mod, control = Rprobit_o$control)
    #  Rprobit_o$H               <- (num_Hessian + t(num_Hessian))/2
    #}


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
