#' calculate spline function
#'
#' @description
#' This function evaluates the natural cubic splines using knots at the points provided.
#'
#' @param grid_points 
#' A real matrix. Each column corresponds to a grid point. 
#' @param knots 
#' A real vector of knots for the spline. 
#'
#' @return
#' A real matrix. Each column corresponds to a separate grid point. 
#'
#' @keywords internal

cal_spline_basis <- function(grid_points,knots){
  L <- dim(knots)[1] # number of layers/dimensions. 
  
  K <- dim(knots)[2] # number of knots
  B <- dim(grid_points)[2] # number of grid_points
  
  if (L==1){ # one dimensional splines
    x_grid <- grid_points[1,]
    spline_basis <- matrix(0,K,B)
    spline_basis[1,] <- 1 #first row: constant
    spline_basis[2,] <- x_grid # second row: X. 
    
    # natural cubic splines 
    dKm1 = (Xplus3(x_grid-knots[K-1])- Xplus3(x_grid-knots[K]))/(knots[K]-knots[K-1])
    for (k in 1:(K-2)){
      dkm1 = (Xplus3(x_grid-knots[k])- Xplus3(x_grid-knots[K]))/(knots[K]-knots[k])
      spline_basis[k+2,] <- dkm1 - dKm1 
    }
    
  }
  
  if (L==2){ # two dimensional splines obtained as cross product of one dimensional case. 

      x_grid <- grid_points[1,]
      spline_basis_x <- matrix(0,K,B)
      spline_basis_x[1,] <- 1 #first row: constant
      spline_basis_x[2,] <- x_grid # second row: X. 
      
      knots_x <- knots[1,]
      
      # natural cubic splines in x direction
      dKm1 = (Xplus3(x_grid-knots_x[K-1])- Xplus3(x_grid-knots_x[K]))/(knots_x[K]-knots_x[K-1])
      for (k in 1:(K-2)){
        dkm1 = (Xplus3(x_grid-knots_x[k])- Xplus3(x_grid-knots_x[K]))/(knots_x[K]-knots_x[k])
        spline_basis_x[k+2,] <- dkm1 - dKm1 
      }
      
      y_grid <- grid_points[2,]
      spline_basis_y <- matrix(0,K,B)
      spline_basis_y[1,] <- 1 #first row: constant
      spline_basis_y[2,] <- y_grid # second row: X. 
      knots_y <- knots[2,]
      # natural cubic splines in y direction 
      dKm1 = (Xplus3(y_grid-knots_y[K-1])- Xplus3(y_grid-knots_y[K]))/(knots_y[K]-knots_y[K-1])
      for (k in 1:(K-2)){
        dkm1 = (Xplus3(y_grid-knots_y[k])- Xplus3(y_grid-knots_y[K]))/(knots_y[K]-knots_y[k])
        spline_basis_y[k+2,] <- dkm1 - dKm1 
      }
      
      # put the two directions together 
      spline_basis <- matrix(0,K*K,B)
      
      for (kr in 1:K){
        for (kc in 1:K){
          spline_basis[(kr-1)*K+kc,] <- spline_basis_x[kr,] * spline_basis_y[kc,] 
        }
      }
  }
    
  
  return (spline_basis)
}
