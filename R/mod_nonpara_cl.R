#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{mod_latclass_cl} object contains the probit model specification for latent class model.
#'
#' @export

mod_nonpara_cl <- R6::R6Class("mod_nonpara_cl",
      inherit = mod_cl,
      # private fields are derived from others
      private = list(
            .tot_params = 0, 
            .num_grid_points = 1
      ),
                      
      active = list(
                        
          #' @field tot_params number of parameters for beta
          tot_params = function(value) {
            if (missing(value)) {
                 private$.tot_params
            } else {
                 stop("`tot_params` is read only", call. = FALSE)
            }
          },
          
          #' @field num_grid_points number of grid points 
          num_grid_points = function(value) {
            if (missing(value)) {
              private$.num_grid_points
            } else {
              stop("`num_grid_points` is read only", call. = FALSE)
            }
          }
          
      ),
                      
      public = list(

          #' @field params parameters at the grid points 
          params = matrix(0,1,1), 
          
          #' @description
          #' Set parameters for grid_points.
          #' @param grid_points A matrix providing the location of the grid points. 
          set_grid_points = function(grid_points) {
            if (is.numeric(grid_points)) {
              self$params <- grid_points
              private$.num_grid_points <- dim(grid_points)[2] 
            } else {
              cat("grid_points must be a numeric matrix.")
            }
          },
          
          #' @description prints the object.
          print = function() {
            cat(sprintf("mod grid_points object (with %d grid_points):\n",self$num_grid_points))
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
