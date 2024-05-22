#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_nonpara_splines_cl} object contains the probit model specification for non-parametric mixing using splines. 
#'
#' @export

mod_nonpara_splines_cl <- R6::R6Class("mod_nonpara_splines_cl",
      inherit = mod_nonpara_cl,
      # private fields are derived from others
      private = list(
            .num_knots = 1
      ),
                      
      active = list(
          
          #' @field num_knots number of grid points 
          num_knots = function(value) {
            if (missing(value)) {
              private$.num_knots
            } else {
              stop("`num_knots` is read only", call. = FALSE)
            }
          }
          
      ),
                      
      public = list(

          #' @field knots defines the location of the knots of the spline
          knots = matrix(0,1,1), 
          
          #' @description
          #' Set location of knots.
          #' @param knots A matrix providing the location of the knots. 
          set_knots = function(knots) {
            if (is.numeric(knots)) {
              self$knots <- knots
              private$.num_knots <- dim(knots)[2] 
            } else {
              cat("knots must be a numeric matrix.")
            }
          },
          
          #' @description prints the object.
          print = function() {
            cat(sprintf("mod splines object (with %d grid_points and %d knots):\n",private$.num_grid_points,private$.num_knots))
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
