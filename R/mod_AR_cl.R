#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_AR_cl} object contains the probit model specification for temporally dependent error processes. 
#' The temporal correlation is introduced via an autoregressive system.
#' In the model class the lag length p is stored, as well as the initial conditions.  
#'
#' @export

mod_AR_cl <- R6::R6Class("mod_AR_cl",
      inherit = mod_cl,
      # private fields are derived from others
      private = list(
            .tot_params = 0
      ),
                      
      active = list(
                        
          #' @field tot_params number of parameters for beta
          tot_params = function(value) {
            if (missing(value)) {
                 private$.tot_params
            } else {
                 stop("`tot_params` is read only", call. = FALSE)
            }
          }
      ),
                      
      public = list(
                        
          #' @field lag_length lag length of the autoregression. 
          lag_length = 1,

          #' @field stationary indicates whether the stationary distribution shall be used to start the state process
          stationary = FALSE,
          
          #' @description
          #' Set dimension of state sequence
          #' @param lag_length An \code{integer}, the dimension of the state sequence
          set_lag_length = function(lag_length) {
             if (is.numeric(lag_length)&(lag_length>=0)) {
                self$lag_length <- lag_length
                vecs <- dim(self$HL)[1]
                s <- round(-0.5+sqrt(0.25+2*vecs))
                private$.tot_params <- private$.lthb + private$.lthO + private$.lthL + lag_length + self$alt-2
              } else {
                cat("lag_length must be a non-negative integer.")
              }
            },
          
          #' @description
          #' Set stationary initial distribution
          #' @param stationary A \code{Boolean}
          set_stationary = function(stationary) {
            if (is.logica(stationary)) {
              self$stationary <- stationary
            } else {
              cat("stationary must be a Boolean variable.")
            }
          },
          
          #' @description prints the object.
          print = function() {
            cat("mod AR object:\n")
            cat("\n")
            cat(sprintf("alt: %d, lthb: %d, lRE: %d, lthO: %d, lthL: %d \n",
                        self$alt, private$.lthb, private$.lRE, private$.lthO, private$.lthL))
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
            if (self$stationary){
              cat("AR process is started using the stationary distribution.")
            } else {
              cat("AR process is started as zero.")
            }
            invisible(self)
          }
    )
)
