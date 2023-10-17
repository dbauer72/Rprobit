#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_latclass_cl} object contains the probit model specification for latent class model.
#'
#' @export

mod_latclass_cl <- R6::R6Class("mod_latclass_cl",
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
                        
          #' @field num_class number of classes
          num_class = 2,

          #' @description
          #' Set number of classes.
          #' @param num_class An \code{integer}, the number of classes.
          set_num_class = function(num_class) {
             if (is.numeric(num_class)&(num_class>1)) {
                self$num_class <- num_class
                private$.tot_params <- num_class*(private$.lthb + private$.lthO + private$.lthL + 1) - 1 
              } else {
                cat("num_class must be an integer larger than 1.")
              }
            },
          
          #' @description prints the object.
          print = function() {
            cat(sprintf("mod latclass object (with %d classes):\n",self$num_class))
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
            invisible(self)
          }
    )
)
