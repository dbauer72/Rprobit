#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_cl} object contains the probit model specification.
#'
#' @export

mod_cl <- R6::R6Class("mod_cl",

  private = list(
    .lthb = 0,
    .lthO = 0,
    .lthL = 0,
    .lRE = 0,
    .fb = matrix(0, 0, 0),
    .Hb = matrix(0, 0, 0),
    .HO = matrix(0, 0, 0),
    .fO = matrix(0, 0, 0),
    .HL = matrix(0, 0, 0),
    .fL = matrix(0, 0, 0)
  ),

  active = list(

    #' @field lthb number of parameters for beta
    lthb = function(value) {
      if (missing(value)) {
        private$.lthb
      } else {
        stop("`lthb` is read only", call. = FALSE)
      }
    },

    #' @field lthO number of parameters for Omega
    lthO = function(value) {
      if (missing(value)) {
        private$.lthO
      } else {
        stop("`lthO` is read only", call. = FALSE)
      }
    },

    #' @field lRE number of random beta coefficients
    lRE = function(value) {
      if (missing(value)) {
        private$.lRE
      } else {
        stop("`lRE` is read only", call. = FALSE)
      }
    },

    #' @field  lthL number of parameters for Sigma
    lthL = function(value) {
      if (missing(value)) {
        private$.lthL
      } else {
        stop("`lthL` is read only", call. = FALSE)
      }
    },

    #' @field fb vector of fixed parameters for beta
    fb = function(value) {
      if (missing(value)) {
        private$.fb
      } else {
        stopifnot(is.matrix(value))
        private$.fb <- value
      }
    },

    #' @field Hb matrix for dependence of beta on parameters
    Hb = function(value) {
      if (missing(value)) {
        private$.Hb
      } else {
        stopifnot(is.matrix(value))
        private$.Hb <- value
        private$.lthb <- dim(value)[2]
      }
    },

    #' @field fO vector of fixed params for Omega
    fO = function(value) {
      if (missing(value)) {
        private$.fO
      } else {
        stopifnot(is.matrix(value))
        private$.fO <- value
      }
    },

    #' @field HO matrix for dependence of Omega parameters
    HO = function(value) {
      if (missing(value)) {
        private$.HO
      } else {
        stopifnot(is.matrix(value))
        private$.HO <- value
        private$.lRE <- -0.5 + sqrt(0.25 + 2 * dim(value)[1])
        private$.lthO <- dim(value)[2]
      }
    },

    #' @field fL vector of fixed params for Sigma
    fL = function(value) {
      if (missing(value)) {
        private$.fL
      } else {
        stopifnot(is.matrix(value))
        private$.fL <- value
      }
    },

    #' @field HL matrix for dependence of Sigma parameters
    HL = function(value) {
      if (missing(value)) {
        private$.HL
      } else {
        stopifnot(is.matrix(value))
        private$.HL <- value
        if (self$ordered == FALSE) {
          self$alt <- round(-0.5 + sqrt(2 * dim(self$HL)[1] + 0.25))
        }
        private$.lthL <- dim(value)[2]
      }
    }
  ),

  public = list(

    #' @field alt number of alternatives
    alt = 2,

    #' @field ordered \code{TRUE} for ordered choices and \code{FALSE} else
    ordered = FALSE,

    #' @description  initialization function
    #' @param Tp number of choices per decider
    #' @param Hb matrix
    #' @param fb vector
    #' @param HO matrix
    #' @param fO vector
    #' @param HL matrix
    #' @param fL vector
    #' @param alt number of alternatives
    #' @param ordered \code{TRUE} for ordered choices and \code{FALSE} else
    #' @param validate \code{TRUE} to validate specification
    #' @param check_identifiability \code{TRUE} to check for identifiability
    initialize = function(
      Hb = matrix(0, 0, 0), fb = matrix(0, 0, 0), HO = matrix(0, 0, 0),
      fO = matrix(0, 0, 0), HL = matrix(0, 0, 0), fL = matrix(0, 0, 0),
      alt = 2, ordered = FALSE, validate = TRUE, check_identifiability = FALSE
    ) {

      stopifnot(is.matrix(Hb))
      stopifnot(is.matrix(fb))
      stopifnot(is.matrix(HO))
      stopifnot(is.matrix(fO))
      stopifnot(is.matrix(HL))
      stopifnot(is.matrix(fL))
      stopifnot(is.numeric(alt))
      stopifnot(isTRUE(ordered) || isFALSE(ordered))
      stopifnot(isTRUE(validate) || isFALSE(validate))
      stopifnot(isTRUE(check_identifiability) || isFALSE(check_identifiability))

      private$.lthb <- dim(Hb)[2]
      private$.lthO <- dim(HO)[2]
      private$.lthL <- dim(HL)[2]

      M <- dim(HO)[1]
      private$.lRE <- round(-0.5 + sqrt(2 * M + 0.25))

      self$Hb <- Hb
      self$fb <- fb
      self$HO <- HO
      self$fO <- fO
      self$HL <- HL
      self$fL <- fL
      self$alt <- alt
      if (ordered == FALSE) {
        self$alt <- round(-0.5 + sqrt(2 * dim(self$HL)[1] + 0.25))
      }
      self$ordered <- ordered

      if (validate) self$validate()
      if (check_identifiability) check_identifiability(self)
    },

    #' @description validates if parameters are conformable
    validate = function() {
      if (dim(self$Hb)[1] != dim(self$fb)[1]) {
        warning(
          "Dimensions of Hb and fb do not match. Resizing fb.",
          call. = FALSE
        )
        self$fb <- matrix(0, dim(self$Hb)[1], 1)
      }
      if (dim(self$HO)[1] != dim(self$fO)[1]) {
        warning(
          "Dimensions of HO and fO do not match. Resizing fO.",
          call. = FALSE
        )
        self$fO <- matrix(0, dim(self$HO)[1], 1)
      }
      if (dim(self$HL)[1] != dim(self$fL)[1]) {
        warning(
          "Dimensions of HL and fL do not match. Resizing fL.",
          call. = FALSE
        )
        self$fL <- matrix(0, dim(self$HL)[1], 1)
      }
    },

    #' @description
    #' Set number of alternatives.
    #' @param alt An \code{integer}, the number of alternatives.
    set_alt = function(alt) {
      if (self$ordered == FALSE) {
        cat("For unordered probit models the number of alternatives is calculated based on Sigma.")
      } else {
        if (is.numeric(alt)) {
          self$alt <- alt
        } else {
          cat("alt must be an integer.")
        }
      }
    },

    #' @description prints the object.
    print = function() {
      cat("mod object:\n")
      cat("\n")
      cat(
        sprintf(
          "alt: %d, lthb: %d, lRE: %d, lthO: %d, lthL: %d \n",
          self$alt, private$.lthb, private$.lRE, private$.lthO, private$.lthL
        )
      )
      cat("\n")
      cat("b (Hb, fb): \n")
      print_est(self$Hb)
      print_est(self$fb)
      cat("\n")
      if (private$.lRE > 0) {
        cat("Omega (HO, fO): \n")
        print_est(self$HO)
        print_est(self$fO)
      } else {
        cat("Omega: 0 \n")
      }
      cat("\n")
      cat("Sigma: \n")
      print_est(self$HL)
      print_est(self$fL)
      if (self$ordered) {
        cat(sprintf("Ordered probit with %d alternatives. \n", self$alt))
      }
      invisible(self)
    }
  )
)
