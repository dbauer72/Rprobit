#' R6 Object Representing the Estimation Specification
#'
#' @description
#' A \code{control_cl} object contains the estimation specification.

control_cl <- R6::R6Class("control_cl",
  public = list(

    #' @field control_nlm
    #' A \code{list} of arguments controlling \code{\link[stats]{nlm}}.
    control_nlm = list(),

    #' @field probit
    #' * If \code{TRUE}, the full probit likelihood is used for estimation.
    #' * If \code{FALSE} the CML is used.
    probit = TRUE,

    #' @field approx_method
    #' A \code{character}, indicating which method to use for the approximation
    #' of the normal CDF. Can be:
    #' * \code{"SJ"} (for Solow–Joe approximation)
    #' * \code{"ME"} (for Mendell–Elston approximation)
    #' * \code{"TVBS"} (for Two-Variate Bivariate Screening approximation)
    approx_method = "SJ",

    #' @field hess
    #' * If \code{FALSE}, the calculation of the analytic Hessian is avoided.
    #' * Else, it is calculated.
    hess = FALSE,

    #' @field pairs_list
    #' A \code{list} for using different weighting schemes in the CML.
    pairs_list = NULL,

    #' @field el
    #' * If \code{TRUE}, the empirical likelihood is used
    #' * Else, it is not.
    el = TRUE,

    #' @field control_weights
    #' A \code{list} of weights for the CML.
    control_weights = NULL,

    #' @field nCores
    #' An \code{integer}, the number of cores to use for parallel computation.
    nCores = 1,

    #' @field normalize
    #' * If \code{TRUE}, a normalization of the regressors is added.
    #' * If \code{FALSE}, it is not.
    normalize = FALSE,

    #' @description
    #' initialization function
    #' @param control_nlm
    #' control list
    #' @param probit
    #' Boolean whether the full probit likelihood or a CML should be used
    #' @param approx_method
    #' which approximation method to use. SJ|ME|TVBS.
    #' @param hess
    #' Boolean; if FALSE calculation of the analytic Hessian is avoided.
    #' @param pairs_list
    #' used for using different weighting schemes in the CML
    #' @param el
    #' Boolean whether the empirical likelihood is employed
    #' @param control_weights
    #' a list of control weights
    #' @param nCores
    #' number of cores to use for parallel calculation.
    #' @param normalize
    #' Boolean add a normalization of the regressors?
    initialize = function(
      control_nlm = list(), probit = TRUE, approx_method = "SJ", hess = FALSE,
      pairs_list = NULL, el = FALSE, control_weights = NULL, nCores = 1,
      normalize = FALSE
    ) {
      stopifnot(is.list(control_nlm), is.logical(probit))
      stopifnot(is.logical(hess))
      stopifnot(is.logical(el))
      stopifnot(is.numeric(nCores), length(nCores) == 1, nCores > 0)
      stopifnot(is.logical(normalize) | is.list(normalize))

      self$control_nlm <- control_nlm
      self$probit <- probit
      self$approx_method <- approx_method
      self$hess <- hess
      self$pairs_list <- pairs_list
      self$el <- el
      self$control_weights <- control_weights
      self$nCores <- nCores
      self$normalize <- normalize
    },

    #' @description
    #' Set control_nlm
    #' @param val New control_nlm
    set_control_nlm = function(val) {
      if (is.list(val)) {
        self$control_nlm <- val
      } else {
        cat("control_nlm must be a list.")
      }
    },

    #' @description
    #' Set probit
    #' @param val New probit
    set_probit = function(val) {
      if (is.logical(val)) {
        self$probit <- val
      } else {
        cat("probit must be boolean.")
      }
    },

    #' @description
    #' Set normalize
    #' @param val New normalize
    set_normalize = function(val) {
      if (is.logical(val) | is.list(val)) {
        self$normalize <- val
      } else {
        cat("normalize must be boolean or list.")
      }
    },

    #' @description
    #' Set approx_method
    #' @param val New  approx_method
    set_approx_method = function(val) {
      if (is.character(val)) {
        self$approx_method <- val
      } else {
        self$approx_method <- "SJ"
      }
    },

    #' @description
    #' Set hess
    #' @param val New hess
    set_hess = function(val) {
      if (is.logical(val)) {
        self$hess <- val
      } else {
        cat("hess must be a boolean.")
      }
    },

    #' @description
    #' Set el
    #' @param val New el
    set_el = function(val) {
      if (is.logical(val)) {
        self$el <- val
      } else {
        cat("el must be a boolean.")
      }
    },

    #' @description
    #' Set nCores
    #' @param val New nCores
    set_nCores = function(val) {
      if (is.numeric(val) && (length(val) == 1) && (val > 0)) {
        self$nCores <- val
      } else {
        cat("nCores must be a scalar number.")
      }
    },

    #' @description
    #' print object
    print = function() {
      cat(sprintf("Controls: number cores: %d, approximation method: %s \n", self$nCores, self$approx_method))
      if (self$probit) {
        cat("Estimated: full probit likelihood.\n")
      } else {
        cat("Estimated: CML.")
      }
      cat(sprintf("Features: hess: %d, empirical likelihood: %d\n", self$hess, self$el))
    }
  )
)
