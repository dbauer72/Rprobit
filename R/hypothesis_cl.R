#' R6 Object Representing hypothesis Specifications
#'
#' @description
#' A \code{hypothesis_cl} object contains the specification of a linear hypothesis for a probit model. 
#'
#' @export

hypothesis_cl <- R6::R6Class("hypothesis_cl",
     public = list(
                               
       #' @field Hb
       #' A real matrix storing information for restriction beta = Hb*theta + fb.
       Hb = NULL, 
       
       #' @field fb
       #' A real matrix storing information for restriction beta = Hb*theta + fb.
       fb = NULL, 
       
       #' @field HO
       #' A real matrix storing information for restriction vOm = HO*theta_O + fO.
       HO = NULL, 
       
       #' @field fO
       #' A real matrix storing information for restriction vOm = HO*theta_O + fO.
       fO = NULL, 
       
       #' @field HL
       #' A real matrix storing information for restriction vL = HL*theta_L + fL.
       HL = NULL, 
       
       #' @field fL
       #' A real matrix storing information for restriction vL = HL*theta_L + fL.
       fL = NULL, 
       
       #' @description  initialization function
       #' @param Hb matrix 
       #' @param fb vector
       #' @param HO matrix
       #' @param fO vector
       #' @param HL matrix
       #' @param fL vector
       initialize = function(Hb = matrix(0, 0, 0), fb = matrix(0, 0, 0), HO = matrix(0, 0, 0),
                fO = matrix(0, 0, 0), HL = matrix(0, 0, 0), fL = matrix(0, 0, 0)) {
           
           stopifnot(is.matrix(Hb))
           stopifnot(is.matrix(fb))
           stopifnot(is.matrix(HO))
           stopifnot(is.matrix(fO))
           stopifnot(is.matrix(HL))
           stopifnot(is.matrix(fL))
           
           self$fb <- fb
           if (dim(fb)[1] == dim(Hb)[1]) {
             self$Hb <- Hb
           } else {
             self$Hb = matrix(0,dim(fb)[1],0)
           }
           self$fO <- fO
           if (dim(fO)[1] == dim(HO)[1]) {
             self$HO <- HO
           } else {
             self$HO = matrix(0,dim(fO)[1],0)
           }
           self$fL <- fL
           if (dim(fL)[1] == dim(HL)[1]) {
             self$HL <- HL
           } else {
             self$HL = matrix(0,dim(fL)[1],0)
           }         
       }, 
       
       #' @description validates if parameters are conformable.
       validate = function() {
         # checks if dimensions match
         if (dim(self$Hb)[1] != dim(self$fb)[1]) {
           warning("Dimensions of Hb and fb do not match. Resizing fb")
           self$fb <- matrix(0, dim(self$Hb)[1], 1)
         }
         if (dim(self$HO)[1] != dim(self$fO)[1]) {
           warning("Dimensions of HO and fO do not match. Resizing fb")
           self$fO <- matrix(0, dim(self$HO)[1], 1)
         }
         if (dim(self$HL)[1] != dim(self$fL)[1]) {
           warning("Dimensions of HL and fL do not match. Resizing fb")
           self$fL <- matrix(0, dim(self$HL)[1], 1)
         }
       }
       
     )                 
)
