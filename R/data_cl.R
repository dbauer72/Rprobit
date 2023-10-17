#' R6 code for setting up a  'data'-object
#' @description \code{data_cl}-object.
#' @param mod \code{data_cl}-object
#' @return (valid) \code{data_cl}-object

data_cl <- R6::R6Class("data_cl",
  private = list(
    .N = 0,
    .Tp = 0
  ),
  active = list(
    #' @field N number of deciders
    N = function(value) {
      if (missing(value)) {
        private$.N
      } else {
        stop("`N` is read only", call. = FALSE)
      }
    },
    #' @field Tp choice occasions per decider
    Tp = function(value) {
      if (missing(value)) {
        private$.Tp
      } else {
        stop("`Tp` is read only", call. = FALSE)
      }
    }
  ),
  public = list(
    #' @field data list containing the data in regressor form
    data = list(),
    #' @field class_member list  membership for latent classes 
    class_member = list(), 
    
    #' @field ordered Boolean indicating if choices are ordered.
    ordered = FALSE,
    #' @field vars string array containing the names of the regressor variables.
    vars = "",
    #' @description  initialization function
    #' @param data list
    #' @param ordered Boolean indicating if choices are ordered.
    #' @param vars strings of regressor names
    initialize = function(data = list(), ordered = FALSE, vars = "") {
      stopifnot(is.list(data), is.character(vars), is.logical(ordered))
      self$data <- data
      N <- length(data)
      private$.N <- N
      private$.Tp <- rep(0, N)
      if (N > 0) {
        for (j in 1:N) {
          private$.Tp[j] <- length(data[[j]]$y)
        }
      }
      self$vars <- vars
      self$ordered <- ordered
    },
    #' @description
    #' Set data
    #' @param val data frame
    set_data = function(val) {
      if (is.list(val)) {
        self$data <- val
        N <- length(val)
        private$.N <- N
        private$.Tp <- rep(0, N)
        for (j in 1:N) {
          private$.Tp[j] <- length(val[[j]]$y)
        }
      } else {
        cat("data must be a list.")
      }
    },
    #' @description
    #' Set class_member
    #' @param val list
    set_class_member = function(val) {
      if (is.list(val)) {
        self$class_member <- val
      } else {
        cat("class_member must be a list.")
      }
    },
    #' @description
    #' Set vars
    #' @param val vars
    set_vars = function(val) {
      if ((is.character(val)) && (length(val) == dim(data[[1]]$X[[1]])[2])) {
        self$vars <- val
      } else {
        cat("'vars' must be a string array assigning one label to each regressor.")
      }
    },
    #' @description
    #' Set ordered
    #' @param val ordered
    set_ordered = function(val) {
      if (is.logical(val)) {
        self$ordered <- val
      } else {
        cat("ordered must be a Boolean.")
      }
    },
    #' @description print object
    print = function() {
      cat(sprintf(" N: %d, Tp: from %d to %d \n", private$.N, min(private$.Tp), max(private$.Tp)))
      cat(sprintf("Variables: %s \n", paste(self$vars, collapse = ", ")))
    },
    #' @description plot the data frame: density plot for each chosen alternative.
    #'  plot the data set
    plot = function() {
      y1 <- self$data[[1]]$y[1]
      X1 <- self$data[[1]]$X[[1]]

      vary <- colnames(X1)

      df <- data.frame(choice = y1, matrix(X1, nrow = 1))
      if (private$.Tp[1] > 1) {
        for (tp in 2:private$.Tp[1]) {
          y1 <- self$data[[1]]$y[tp]
          X1 <- self$data[[1]]$X[[tp]]

          df <- rbind(df, data.frame(choice = y1, matrix(X1, nrow = 1)))
        }
      }

      if (private$.N > 1) {
        for (n in 2:private$.N) {
          for (tp in 1:private$.Tp[n]) {
            y1 <- self$data[[n]]$y[tp]
            X1 <- self$data[[n]]$X[[tp]]

            df <- rbind(df, data.frame(choice = y1, matrix(X1, nrow = 1)))
          }
        }
      }
      choice <- df[, "choice"]
      if (is.factor(choice)) {
        alt <- levels(choice)
      } else {
        alt <- unique(choice)
      }

      vars <- "choice"

      for (jv in 1:length(vary)) {
        for (ja in 1:length(alt)) {
          vars[1 + ja + (length(alt)) * (jv - 1)] <- paste(vary[jv], "_", alt[ja], sep = "")
        }
      }

      colnames(df) <- vars

      bars <- rep(0, length(vary))
      tab <- list()
      for (jv in 1:length(vary)) {
        xv <- c()
        for (ja in 1:length(alt)) {
          xv <- c(xv, df[, paste(vary[jv], "_", alt[ja], sep = "")])
        }
        if (length(unique(xv)) < 11) {
          bars[jv] <- 1
        }
      }

      par(mfrow = c(length(alt), length(vary)))
      for (ja in 1:length(alt)) {
        ind <- which(choice == alt[ja])
        for (jv in 1:length(vary)) {
          if (bars[jv] == 0) {
            xv <- df[ind, paste(vary[jv], "_", alt[2], sep = "")]
            plot(density(xv), main = vary[jv], xlab = paste("Choice: ", alt[ja]), col = 2)
            for (jc in 3:length(alt)) {
              xv <- df[ind, paste(vary[jv], "_", alt[jc], sep = "")]
              if (var(xv) > 0) {
                lines(density(xv), col = jc)
              }
            }
          } else {
            xv <- data.frame(x = df[ind, paste(vary[jv], "_", alt[2], sep = "")], a = alt[2])
            for (jc in 3:length(alt)) {
              xv <- rbind(xv, data.frame(x = df[ind, paste(vary[jv], "_", alt[jc], sep = "")], a = alt[jc]))
            }
            tt <- table(xv$x, xv$a)

            for (t in 1:dim(tt)[2]) {
              tt[, t] <- tt[, t] / sum(tt[, t])
            }
            barplot(t(tt), beside = TRUE, col = 2:length(alt), main = vary[jv], xlab = paste("Choice: ", alt[ja]))
          }
        }
      }

      par(mfrow = c(1, 1))
    }
  )
)
