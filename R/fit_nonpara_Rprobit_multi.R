#' Fit a Probit Model using non-parametric mixing distribution based on a fixed grid.
#'
#' @description
#' This function fits a mixed probit model using a non-parametric mixing with fixed grid and estimated mixing
#' probabilities.
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
#' @param cml_pair_type
#' An integer.
#' @param iter.adapt
#' an integer specifying the number of adaptations. iter.adapt (default) implies none.
#' cml_pair_type = 0 (full pairwise), = 1 (adjacent pairwise looped), = -1 (single choices)
#' By default, \code{cml_pair_type = 1}.
#'
#'
#' @return
#' A \code{\link{Rprobit_cl}} object.
#'
#' @export

fit_nonpara_Rprobit_multi <- function(Rprobit_obj, init_method = "random", control_nlm = NULL,cml_pair_type  = 1, iter.adapt = 1) {

  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  mod_orig <- Rprobit_o$mod$clone(deep=TRUE)
  data_raw_orig = Rprobit_o$data_raw$clone(deep=TRUE)

  #####################################################
  ## set up helper functions for optimization.    #####
  #####################################################
  # helper functions needed in the algorithm


  ### extract print level
  print.level           <- 0

  if(!is.null(control_nlm$print.level)){
    print.level    <- control_nlm$print.level
  } else if(!is.null(Rprobit_o$control$control_nlm$print.level)){
    print.level    <- Rprobit_o$control$control_nlm$print.level
  }




  ### define optimizer function
  optimizer <- function(ll_function, theta, Rprobit_o, data, cml_pair_type){


    ### define nlm parameters
    p_nlm <- theta


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
    browser()

    ### optimise
    out <- stats::nlm(f                  = ll_function,
                      p                  = p_nlm,
                      data               = data,
                      mod                = Rprobit_o$mod,
                      control            = Rprobit_o$control,
                      cml_pair_type      = cml_pair_type,
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




  ### define the likelihood function
  ll_function <- function(theta, data_tr, mod,control, cml_pair_type) {

    mod_new <- mod$clone(deep = TRUE)

    #browser()

    params <- matrix(0, length(theta) + 1, mod_new$num_grid_points)
    params[1,] <- mod_new$params
    for (j in 1:mod_new$num_grid_points){
      params[-1,j] <- theta
    }
    mod_new$params <- params
    
    out <- fit_nonpara_grid(data_tr, mod_new, control, cml_pair_type  = cml_pair_type)
    cr <- - (out$weights %*% log(out$lprob))
    
    return(cr)
  }


  ### Check if parallel framework should be used and transform data in that case
  use_parallel <- FALSE
  if (!is.null(Rprobit_o$control$nCores)) {
    if (Rprobit_o$control$nCores > 1) {
      use_parallel <- TRUE
    }
  }

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
  K <- Rprobit_o$mod$num_grid_points
  tot_params <- K

  ### perform numerical optimization and save results
  Rprobit_o$control$hess <- FALSE
  if (print.level > 0) {
    print("Start optimisation...")
  }

  Rprobit_o$data = data_tr

  num_params <- Rprobit_o$mod$lthb + Rprobit_o$mod$lthO + Rprobit_o$mod$lthL - 1
  theta = rnorm(num_params)

  time_1 <- Sys.time()
  nlm_out         <- optimizer(ll_function = ll_function, theta, Rprobit_o = Rprobit_o, data = Rprobit_o$data, cml_pair_type)

  time_2          <- Sys.time()

  #Rprobit_o$mod$params <- grid_points
  Rprobit_o$theta <- nlm_out$estimate

  neg_ll_fit <- ll_function(Rprobit_o$theta, Rprobit_o$mod,Rprobit_o$control, cml_pair_type)
  Rprobit_o$ll <- (-1) * neg_ll_fit

  Rprobit_o$fit <- approx_method
  Rprobit_o$control$hess <- hessian
  if (!is.null(Rprobit_o$info)) {
    Rprobit_o$info$estimation_time <- difftime(time_2, time_1, units = "secs")
    Rprobit_o$info$nlm_info <- nlm_out
  } else {
    Rprobit_o$info <- list(
      estimation_time = difftime(time_2, time_1, units = "secs"),
      nlm_info <- nlm_out
    )
  }

  ### computation of variance


  time_1                    <- Sys.time()

  Rprobit_o$H               <- diag(length(theta)) #out_opt$hessian
  Rprobit_o$J               <- diag(length(theta)) # out_opt$hessian

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
  #Rprobit_o$mod$params <- grid_points
  Rprobit_o$data_raw <- data_raw_orig

  return(Rprobit_o)
}
