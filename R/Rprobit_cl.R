#' R6 Object Representing a Probit Model
#'
#' @description
#' An \code{Rprobit_cl} object contains the probit model specification.

Rprobit_cl <- R6::R6Class("Rprobit_cl",
  public = list(

    #' @field control
    #' An \code{\link{control_cl}} object which stores information for estimation.
    control = NULL,

    #' @field data
    #' A \code{\link{data_cl}} object.
    data = NULL,

    #' @field data_raw
    #' A \code{\link{data_raw_cl}} object.
    data_raw = NULL,

    #' @field form
    #' A \code{formula} object.
    form = NULL,

    #' @field re
    #' A \code{character} vector containing the variable names for random effects.
    re = NULL,

    #' @field mod
    #' A \code{\link{mod_cl}} object.
    mod = NULL,

    #' @field vars
    #' A \code{character} vector containing variable names.
    vars = NULL,

    #' @field alt_names
    #' A \code{character} vector containing alternative names.
    alt_names = NULL,

    #' @field theta
    #' A \code{numeric} parameter vector.
    theta = NULL,

    #' @field theta_0
    #' A \code{numeric} true parameter vector (for simulated data only).
    theta_0 = NULL,

    #' @field ll
    #' A \code{numeric}, the log likelihood value.
    ll = NA,

    #' @field H
    #' A \code{matrix}, the Hessian at estimate.
    H = NA,

    #' @field J
    #' A variance \code{matrix} of score vector.
    J = NA,

    #' @field grad
    #' A gradient \code{vector} at estimate.
    grad = NA,

    #' @field fit
    #' A \code{character}, the approximation method.
    fit = NULL,

    #' @field vv
    #' A variance \code{matrix} computed using the Hessian.
    vv = NULL,

    #' @field vv2
    #' A variance \code{matrix} of different variance estimate computed using
    #' the approximated Hessian.
    vv2 = NULL,

    #' @field info
    #' A \code{list} containing information on fitting process.
    info = NULL,
    #' @field gradEL list containing gradient contribution of each individual
    gradEL = NULL,
    #' @field HessEL list containing Hessian contribution of each individual
    HessEL = NULL,
    
    #' @description  initialization function
    #' @param data      object of class \link{data_cl}
    #' @param data_raw  object of class \link{data_raw_cl}
    #' @param form      \link{formula} object
    #' @param re        list of strings containing the variable names for random effects
    #' @param mod       object of class \link{mod_cl}
    #' @param vars      list of strings containing variable names
    #' @param alt_names list of strings of names for alternatives
    #' @param theta     parameter vector
    #' @param theta_0   true parameter vector for simulated data
    #' @param ll        log likelihood value
    #' @param H         matrix of Hessian at parameter estimate
    #' @param J         matrix of variance of score vector
    #' @param grad      gradient vector at estimate
    #' @param fit       approximation method
    #' @param vv        matrix of variance computed using the Hessian
    #' @param vv2       matrix of different variance estimate computed using the approximated Hessian
    #' @param control   object of class \link{control_cl}. Stores all information for estimation.
    #' @param info      list containing information on fitting process
    #' @param gradEL    list containing gradient contribution of each individual
    #' @param HessEL    list containing Hessian contribution of each individual
    initialize = function(data = NULL,
                          data_raw = NULL,
                          form = NULL,
                          re = NULL,
                          mod = NULL,
                          vars = NULL,
                          alt_names = NULL,
                          theta = NULL,
                          theta_0 = NULL,
                          ll = NA,
                          H = NA,
                          J = NA,
                          grad = NA,
                          fit = NA,
                          vv = NA,
                          vv2 = NA,
                          control = NULL,
                          info = NULL,
                          gradEL = NULL,
                          HessEL = NULL) {
      if (inherits(data, "data_cl")) {
        self$data <- data
      } else {
        self$data <- NULL
      }
      if (inherits(data_raw, "data_raw_cl")) {
        self$data_raw <- data_raw
      } else {
        self$data_raw <- NULL
      }
      if (inherits(form, "formula")) {
        self$form <- form
      } else {
        self$form <- NULL
      }
      self$re <- re
      if (inherits(mod, "mod_cl")||inherits(mod, "mod_latclass_cl")) {
        self$mod <- mod
      } else {
        self$mod <- NULL
      }
      self$vars <- vars
      self$alt_names <- alt_names
      self$theta <- theta
      self$theta_0 <- theta_0
      self$ll <- ll
      self$H <- H
      self$J <- J
      self$grad <- grad
      self$fit <- fit
      self$vv <- vv
      self$vv2 <- vv2
      if (inherits(control, "control_cl")) {
        self$control <- control
      } else {
        self$control <- NULL
      }
      self$info <- info
      self$gradEL <- gradEL
      self$HessEL <- HessEL
    },

    #' @description
    #' Set data
    #' @param val New data
    set_data = function(val) {
      if (inherits(val, "data_cl")) {
        self$data <- val
      } else {
        cat("data must be a data_cl object.")
      }
    },

    #' @description
    #' Set data_raw
    #' @param val New data_raw
    set_data_raw = function(val) {
      if (inherits(val, "data_raw_cl")) {
        self$data_raw <- val
      } else {
        cat("data_raw must be a data_raw_cl object.")
      }
    },

    #' @description
    #' Set form
    #' @param val form
    set_form = function(val) {
      if (inherits(val, "formula")) {
        self$form <- val
      } else {
        cat("form must be a formula object.")
      }
    },

    #' @description
    #' Set re
    #' @param val list of strings containing names of random coefficients
    set_re = function(val) {
      if (is.character(val) == TRUE) {
        self$re <- val
      } else {
        cat("re must be characters.")
      }
    },

    #' @description
    #' Set mod
    #' @param val object of class 'mod_cl'
    set_mod = function(val) {
      if (inherits(val, "mod_cl")||inherits(mod, "mod_lat_class_cl")) {
        self$mod <- val
      } else {
        cat("mod must be a mod_cl or mod_latclass_cl object.")
      }
    },

    #' @description
    #' Set vars
    #' @param val list of strings containing names of variables
    set_vars = function(val) {
      if (is.character(val) == TRUE) {
        self$vars <- val
      } else {
        cat("vars must be characters.")
      }
    },

    #' @description
    #' Set alt_names
    #' @param val list of strings containing names of alternatives
    set_alt_names = function(val) {
      if (is.character(val) == TRUE) {
        self$alt_names <- val
      } else {
        cat("alt_names must be characters.")
      }
    },

    #' @description
    #' Set theta
    #' @param val numeric vector
    set_theta = function(val) {
      if (is.numeric(val) == TRUE) {
        self$theta <- val
      } else {
        cat("theta must be characters.")
      }
    },

    #' @description
    #' Set theta_0
    #' @param val numeric vector
    set_theta_0 = function(val) {
      if (is.numeric(val) == TRUE) {
        self$theta_0 <- val
      } else {
        cat("theta_0 must be characters.")
      }
    },

    #' @description
    #' Set ll
    #' @param val double
    set_ll = function(val) {
      if ((is.numeric(val) == TRUE) & (length(ll) == 1)) {
        self$ll <- val
      } else {
        cat("ll must be real")
      }
    },

    #' @description
    #' Set H
    #' @param val Hessian matrix
    set_H = function(val) {
      if (is.matrix(val) == TRUE) {
        self$H <- val
      } else {
        cat("H must be a real marix.")
      }
    },

    #' @description
    #' Set J
    #' @param val varianc matrix of score
    set_J = function(val) {
      if (is.matrix(val) == TRUE) {
        self$J <- val
      } else {
        cat("J must be a real marix.")
      }
    },


    #' @description
    #' Set grad
    #' @param val numeric vector
    set_grad = function(val) {
      if (is.vector(val) == TRUE) {
        self$grad <- val
      } else {
        cat("grad must be a vector.")
      }
    },

    #' @description
    #' Set fit
    #' @param val fit
    set_fit = function(val) {
      # if (is.numeric(val)==TRUE){
      self$fit <- val
      # } else {
      #  cat("fit must be ??")
      # }
    },

    #' @description
    #' Set vv
    #' @param val variance matrix
    set_vv = function(val) {
      if (is.matrix(val) == TRUE) {
        self$vv <- val
      } else {
        cat("vv must be a real marix.")
      }
    },

    #' @description
    #' Set vv2
    #' @param val approximate variance matrix
    set_vv2 = function(val) {
      if (is.matrix(val) == TRUE) {
        self$vv2 <- val
      } else {
        cat("vv2 must be a real marix.")
      }
    },


    #' @description
    #' Set control
    #' @param val control object of class \link{control_cl}
    set_control = function(val) {
      if (inherits(val, "control_cl")) {
        self$control <- val
      } else {
        cat("control must be a 'control_cl' object.")
      }
    },

    #' @description
    #' Set info
    #' @param val info
    set_info = function(val) {
      # if (is.numeric(val)==TRUE){
      self$info <- val
      # } else {
      #  cat("fit must be ??")
      # }
    },
    
    #' @description
    #' Set gradEL
    #' @param val list
    set_gradEL = function(val) {
      #if (is.list(val) == TRUE) {
      self$gradEL <- val
      #} else {
      #  cat("gradEL must be a List.")
      #}
    },

    #' @description
    #' Set HessEL
    #' @param val list
    set_HessEL = function(val) {
      #if (is.list(val) == TRUE) {
      self$HessEL <- val
      #} else {
      #  cat("HessEL must be a List.")
      #}
    },

    #' @description
    #' print object
    print = function() {
      cat("Rprobit object\n")
      invisible(self)
    },

    #' @description
    #' return coefficients for beta parameters
    #' @return beta
    coef = function() {
      if (!is.null(self$data)) {
        labels <- self$data$vars
      } else {
        labels <- NULL
      }
      par <- build_par_from_mod(self$theta, self$mod, variances = NA, labels = labels)
      return(par$b)
    },

    #' @description
    #' provide plots for model (useful for latent class models)
    #' @param dims 
    #' vector of integers indicating the coordinates to be plotted. Only the first two integers are considered
    #' @param cdf 
    #' Boolean, indicating whether the cdf (TRUE) or the pdf (FALSE; default) should be plotted.  
    #' @param margin
    #' real; factor (0.1 = 10%; minimum default margin) to enlarge plotting area. 
    plot = function(dims = 1, cdf = FALSE, margin = 0.1) {
      plot_Rprobit(self, dims = dims, cdf = cdf, margin = margin)
    },

    #' @description
    #' provide predictions for model contained in object.
    #' @param newdata 
    #' either a \code{\link{data_raw_cl}} object or a data frame. 
    predict = function(newdata = NULL) {
      data_new = self$data_raw$clone()
      if (is.data.frame(newdata)){
          data_new$set_df(newdata)
        } 
      pr <- predict_Rprobit(self, data_new = data_new, all_pred = TRUE)
      return(pr)
    },

    #' @description
    #' provide residuals for model contained in object.
    residuals = function() {
      pr <- predict_Rprobit(self, data_new = self$data_raw, all_pred = TRUE)
      cols <- dim(pr)[2]
      choice <- pr[, cols - 1]
      y <- pr[, 1:(cols - 2)] * 0

      for (nt in 1:length(choice)) {
        y[nt, choice[nt]] <- 1
      }

      res <- y - pr[, 1:(cols - 2)]

      return(res)
    },

    #' @description
    #' provide AIC value for model contained in object.
    AIC = function() {
      #if (is.null(self$theta)) {
      #  out <- NA
      #} else {
        out <- -2 * self$ll + 2 * length(self$theta)
      #}
      return(out)
    },

    #' @description
    #' provide BIC value for model contained in object. Sample size equal to N*T in balanced case
    BIC = function() {
      #if (is.null(self$theta)) {
      #  out <- NA
      #} else {
        out <- -2 * self$ll + log(sum(self$data_raw$Tp)) * length(self$theta)
      #}
      return(out)
    },

    #' @description
    #' provide confusion matrix for model contained in object.
    conf_mat = function() {
      out <- confmat_Rprobit(self)
      return(out$conf_mat)
    }
  )
)
