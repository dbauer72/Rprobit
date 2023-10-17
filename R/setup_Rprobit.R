#' Create Rprobit-object
#'
#' @description
#' This function creates an \code{\link{Rprobit_cl}} object.
#'
#' @inheritParams read_formula
#' @param data_raw
#' \code{data_raw}-object
#' @param mod
#' A \code{\link{mod_cl}} or \code{\link{mod_latclass_cl}} object.
#' @param norm_alt
#' An \code{integer} that labels the alternative which's utility is set to zero
#' (for normalization).
#' @param ids
#' A \code{character} vector, specifying the columns that identify a decision
#' maker uniquely.
#' @param re
#' A \code{character} vector of variable names which are considered to have
#' random effects.
#' @param corr_re
#' A \code{logical}, if \code{TRUE} a full covariance matrix Omega will be
#' setup, if \code{FALSE} (default) only a diagonal one.
#' @param error_struc
#' A \code{character}, specifying the error structure, one of:
#' * \code{"diag"} diagonal matrix Sigma
#' * \code{"fixed"} fixed matrix Sigma
#' * \code{"ident"} Sigma proportional to identity matrix.
#' * \code{"full"} (default)
#' @param seed
#' Seed for the simulation of choice data.
#' No seed set per default.
#' @param theta_0
#' A \code{numeric} vector, the true parameter vector for simulated data.
#' It is drawn randomly if not supplied.
#' @param name
#' A \code{character}, a name for the model.
#' @param description
#' A \code{character}, a description for the model.
#' @param control
#' A \code{list}, contains control variables for estimation procedure.
#' @param avail
#' A \code{list}, contains information on the availability for each alternative
#' and each choice occasion.
#' @param info
#' A \code{list}, contains information on the model and estimation performance.
#' @param alt_names
#' A \code{character} vector, containing the labels for the choices.
#'
#' @return
#' An \code{\link{Rprobit_cl}} object.
#'
#' @export

setup_Rprobit <- function(
    form, data_raw = NULL, mod = NULL, norm_alt = 1, ids = NULL, re = NULL,
    corr_re = FALSE, error_struc = "full", seed = NULL, theta_0 = NULL,
    name = NULL, description = NULL, control = NULL, avail = NULL, info = NULL,
    alt_names = NULL
  ){

  ### save inputs for reproducible
  setup_input <- list(
    form = form, norm_alt = norm_alt, ids = ids, re = re, corr_re = corr_re,
    error_struc = error_struc, seed = seed, theta_0 = theta_0
  )

  if (inherits(control, "control_cl")) {
    control_obj = control
  } else  {
    control_obj = control_cl$new()
    control_obj$control_simulation <- control
  }

  ### check if data_raw or mod are specified
  if (is.null(data_raw) & is.null(mod)) {
    stop("Either 'data_raw' or 'mod' has to be specified.", call. = FALSE)
  }
  if (is.null(data_raw)) {
    if (is.null(mod$alt)) {
      stop("Number 'alt' of choice alternatives has to specified in 'mod'.",
           call. = FALSE)
    }
  }

  ### read 'form'
  read_formula_out  <- read_formula(form)
  ASC               <- read_formula_out$ASC
  allvars           <- read_formula_out$allvars
  choice            <- all.vars(form)[1]

  ### check 'error_struc' and 're'
  if (!(error_struc %in% c("full", "diag", "fixed", "ident"))) {
    warning("Only valid specifications for the error terms are 'diag', 'fixed', 'ident' and 'full'. Proceding with correlated error terms ('full').")
    error_struc <- "full"
  }
  if(!all(re %in% c("ASC",unlist(allvars)))) stop("The following names supplied to 're' are no columns in 'data_raw': \n", paste(re[!(re %in% c("ASC",unlist(allvars)))],collapse = ", "))

  ### extract information from empirical data
  if(!is.null(data_raw)){
    if (inherits(data_raw,"data_raw_cl")){
      data_raw_obj = data_raw
      ids = data_raw$id
      # data_raw already an object
    } else {
      ### create unique identifier
      if(is.null(ids)) {
        stop("'ids' has to be specified.")
      }

      if (is.null(alt_names)){
        alt_names = levels(as.factor(data_raw[,choice]))
      }
      data_raw_obj = data_raw_cl$new(
        df = data_raw,
        alt_names = alt_names,
        id = ids,
        choice = choice,
        varying = allvars[[1]],
        dec_char = allvars[[2]]
      )

      ### save submitted data_raw
      data_raw_0 = data_raw
    }
    ### convert to data from data_raw
    data <- data_raw_to_data(data_raw = data_raw_obj, allvars = allvars, choice = choice, re = re, norm_alt = norm_alt, alt = mod$alt, ASC = ASC)
    
    ### check 'data_raw' and compute number of alternatives
    if (is.null(mod)) {
      mod <- mod_cl$new()

      ualt = unique(data_raw_obj$df[,"choice"])
      alt = length(ualt)

      # parameters for beta
      ltype1 = 0
      if (length(allvars[[1]])==1){
        if (allvars[[1]]==0){
          ltype1 = 0 
        } else {
          ltype1 = 1 
        } 
      } else {
        ltype1 = length(allvars[[1]])
      }
      
      ltype2 = length(allvars[[2]])
      if (allvars[[3]] == 0){
        ltype3=0
        } else {
          ltype3 = length(allvars[[3]])
        }
      lthb <- ltype1  + (ltype2 + ASC) * (alt-1) + ltype3 * alt
      Hb = diag(lthb)
      fb = matrix(0,lthb,1)
      mod$Hb = Hb
      mod$fb = fb


      # random parameters
      lRE <- sum(re %in% allvars[[1]]) + (sum(re %in% allvars[[2]]) + "ASC" %in% re) * (alt - 1) + sum(re %in% allvars[[3]]) * alt
      MO <- lRE * (lRE + 1) / 2
      HO <- diag(MO)
      if (lRE > 0) {
        if (corr_re == TRUE) {
          HO <- diag(MO)
          fO <- matrix(0, MO, 1)
        } else {
          HO <- matrix(0, MO, lRE)
          cur <- 1
          for (j in 1:lRE) {
            HO[cur, j] <- 1
            cur <- cur + lRE - j + 1 ### TODO: try fixing this matrix build
          }
        }
        fO <- matrix(0, MO, 1)
        mod$HO <- HO
        mod$fO <- fO
      }

      # set up error structure
      nL = alt
      HL <- diag((nL*nL-nL)/2+nL)
      if(error_struc == "full"){
        HL    <- HL[, -c(1:(nL+1)), drop=FALSE]
        fL    <- rep(0,(nL*nL-nL)/2+nL)
        fL[nL+1] <- 1

        #lthL  <- sum(HL)
      }
      if(error_struc == "diag"){
        HL    <- HL[, fdiag(nL)+1, drop=FALSE]
        HL    <- HL[, -1, drop=FALSE]
        fL    <- rep(0,nrow(HL))
        #fL    <- rep(0,ncol(HL))
        fL[1] <- 1
        #lthL  <- sum(HL)
      }
      if(error_struc == "ident"){
        HL    <- as.matrix(apply(HL[, fdiag(nL)+1, drop=FALSE],1,sum))
        #HL    <- HL[, -1, drop=FALSE]
        fL    <- rep(0,nrow(HL))
        #fL    <- rep(0,ncol(HL))
        #fL[1] <- 1
        #lthL  <- 1
      }  
        
      if(error_struc == "fixed"){ 
        HL    <- HL[, fdiag(nL)+1, drop=FALSE]
        fL    <- HL %*% rep(sqrt(0.5),ncol(HL))
        HL    <- matrix(0,(nL*nL-nL)/2+nL,0)
        #lthL  <- 0
      }
      mod$HL <- HL
      mod$fL <- as.matrix(fL, ncol = 1)
      # $lthL <- lthL
    }
  }

  ### simulate data_raw
  if (is.null(data_raw)) {
    sim_data_raw_out <- sim_data_raw(
      ASC = ASC,
      allvars = allvars,
      re = re,
      choice = choice,
      theta = theta_0,
      mod = mod,
      seed = seed,
      control_simulation = control_obj$control_simulation
    )
    data_raw <- sim_data_raw_out$data_raw
    if (is.null(alt_names)) {
      alt_names <- levels(as.factor(data_raw$df[, data_raw$choice]))
    }
    data_raw_obj = data_raw_cl$new(df = data_raw$df,alt_names = alt_names,id = "id_macml",choice = choice,varying =allvars[[1]], dec_char = allvars[[2]])


    theta_0 <- sim_data_raw_out$theta_0

    ### save simulated data_raw
    data_raw_0 <- data_raw
    
    ### save simulated data
    data <-  sim_data_raw_out$data
  }


  ### Normalise regressors in data_raw

  if (mod$ordered == FALSE) {
    data_raw_obj <- normalise_data(data_raw = data_raw_obj, normalisation = control_obj$normalize, allvars = allvars, ASC = ASC)
    control_obj$normalize <- attr(x = data_raw_obj, which = "normalisation")


    ### create 'data' from 'data_raw'
    # data = data_raw_to_data(data_raw = data_raw_obj, allvars = allvars, choice = choice, re = re, norm_alt = norm_alt, alt = mod$alt, ASC = ASC)
  } else {
    # for ordered data there is no difference between data_raw and data
    data = data_cl$new(ordered = TRUE,vars = data_raw_obj$dec_char)
    ids = data_raw_obj$df[,data_raw_obj$id]
    unique_ids = unique(ids)
    data_df = list()
    for (j in 1:length(unique_ids)){
      ind = which(ids == unique_ids[j])
      data_df[[j]] = list(X = data.matrix(data_raw_obj$df[ind,data_raw_obj$dec_char]), y = data.matrix(data_raw_obj$df[ind,data_raw_obj$choice]))
    }
    data$set_data(data_df)
  }

  attr(x = data_raw, which = "ids") <- ids

  ### create 'data' from 'data_raw'
  vars <- data$vars

  ### build 'Rprobit_obj'
  Rprobit_obj <- Rprobit_cl$new(
    data      = data,
    data_raw  = data_raw_obj,
    form      = form,
    re        = re,
    mod       = mod,
    alt_names = alt_names,
    vars      = vars,
    theta_0   = theta_0
  )

  if (mod$ordered == FALSE) {
    Rprobit_obj <- check_avail(data_raw = data_raw_obj, avail = avail, Rprobit_obj = Rprobit_obj)
  }

  ### TODO: This line does not work when control already exists!
  Rprobit_obj$control <- control_obj$clone()

  ### add 'info'
  Rprobit_obj$info <- build_info(
    Rprobit_obj = Rprobit_obj,
    data_raw    = data_raw,
    form        = form,
    ids         = ids,
    name        = name,
    description = description,
    control     = control_obj,
    info        = info
  )

  Rprobit_obj$info$setup_input <- setup_input

  return(Rprobit_obj)
}
