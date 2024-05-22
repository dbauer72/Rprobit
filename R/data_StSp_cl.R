#' R6 Object Representing data 
#'
#' @description
#' A \code{data_StSp_cl} object contains a data set in our specific format, extending data_cl.
#' In addition to the data_cl a field 'time' is introduced that lists the times at which samples are taken. 
#'
#' @export

data_StSp_cl <- R6::R6Class("data_StSp_cl",
                                inherit = data_cl,
                                
                                active = list(),
                                
                                public = list(
                                  
                                  #' @field time listing the observation times 
                                  time = c(),
                                  
                                  #' @description
                                  #' Set time vector.
                                  #' @param time A vector, listing the time points
                                  set_time = function(time) {
                                    if (is.numeric(time))  {
                                      self$time <- time
                                    } else {
                                      cat("time must be a numeric vector.")
                                    }
                                  },
                                  
                                  #' @field quest integer indicating the maximal number of questions per person. 
                                  quest = 1, 

                                  #' @description
                                  #' Set question number.
                                  #' @param qu integer; maximal number of questions. 
                                  set_quest = function(qu) {
                                    if (is.numeric(qu))  {
                                      self$quest <- qu
                                    } else {
                                      cat("Quest must be an integer.")
                                    }
                                  },
                                  
                                  
                                  #' @description prints the object.
                                  print = function() {
                                    cat("Data object including time points.")
                                    cat("\n")
                                    
                                    invisible(self)
                                  }
                                )
)
