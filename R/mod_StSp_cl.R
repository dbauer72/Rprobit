#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_StSp_cl} object contains the probit model specification for temporally dependent error processes. 
#' The temporal correlation is introduced via a state space system.
#' In the model class the state dimension n are stored.  
#'
#' @export

mod_StSp_cl <- R6::R6Class("mod_StSp_cl",
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
                        
          #' @field dim_state dimension of state sequence. 
          dim_state = 2,

          #' @field stationary indicates whether the stationary distribution shall be used to start the state process
          stationary = FALSE,
          
          #' @description
          #' Set dimension of state sequence
          #' @param dim_state An \code{integer}, the dimension of the state sequence
          set_dim_state = function(dim_state) {
             if (is.numeric(dim_state)&(dim_state>=0)) {
                self$dim_state <- dim_state
                vecs <- dim(self$HL)[1]
                s <- round(-0.5+sqrt(0.25+2*vecs))
                private$.tot_params <- private$.lthb + private$.lthO + private$.lthL + 2*dim_state*s + self$alt-2
              } else {
                cat("dim_state must be a non-negative integer.")
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
            cat("mod StSp object:\n")
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
              cat("State is started using the stationary distribution.")
            } else {
              cat("State is started as zero.")
            }
            invisible(self)
          }
    )
)
