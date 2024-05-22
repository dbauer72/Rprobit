#' R6 Object Representing Model Specifications
#'
#' @description
#' A \code{data_raw_StSp_cl} object contains a raw data set in the wide format, extending data_raw_cl.
#' In addition to the data_raw_cl a field 'time' is introduced that lists the times at which samples are taken. 
#'
#' @export

data_raw_StSp_cl <- R6::R6Class("data_raw_StSp_cl",
                                inherit = data_raw_cl,
                                # private fields are derived from others
                                private = list(
                                  .time = c()
                                ),
                                
                                active = list(
                                  
                                  #' @field time observation time points
                                  time = function(value) {
                                    if (missing(value)) {
                                      private$.time
                                    } else {
                                      stop("`time` is read only", call. = FALSE)
                                    }
                                  }
                                ),
                                
                                public = list(
                                  
                                  #' @field time_col column of data set containing time points 
                                  time_col = "",
                                  
                                  #' @field quest_col column of data set containing number of question
                                  quest_col = "",
                                  
                                  #' @description
                                  #' Set time column.
                                  #' @param time_col An \code{string}, denoting the column listing the time points
                                  set_time_col = function(time_col) {
                                    if (is.character(time_col)& base::any(base::grepl(time_col,colnames(self$df))))  {
                                      self$time_col <- time_col
                                      private$.time <- base::unique(self$df[time_col]) 
                                    } else {
                                      cat("time_col must be a valid column name.")
                                    }
                                  },
                                  
                                  #' @description
                                  #' Set question column.
                                  #' @param quest_col An \code{string}, denoting the column listing the time points
                                  set_quest_col = function(quest_col) {
                                    if (is.character(quest_col)& any(base::grepl(quest_col,colnames(self$df))))  {
                                      self$quest_col <- quest_col
                                      private$.Tp <- base::max(self$df[quest_col]) 
                                    } else {
                                      cat("quest_col must be a valid column name.")
                                    }
                                  },
                                  
                                  #' @description prints the object.
                                  print = function() {
                                    cat("Raw data object including levels for time points and item numbers.")
                                    cat("\n")
                                    
                                    invisible(self)
                                  }
                                )
)
