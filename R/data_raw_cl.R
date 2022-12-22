#' R6 Object Representing a Data Object
#'
#' @description
#' An \code{data_raw_cl} object contains the data for modeling.

data_raw_cl <- R6::R6Class("data_raw_cl",

  private = list(
    .N = 0,
    .Tp = 0,
    .dec_char = "",
    .varying = ""
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

    #' @field Tp number of choice decisions per decider
    Tp = function(value) {
      if (missing(value)) {
        private$.Tp
      } else {
        stop("`Tp` is read only", call. = FALSE)
      }
    },

    #' @field dec_char decider specific regressors
    dec_char = function(value) {
      if (missing(value)) {
        private$.dec_char
      } else {
        stop("`decider` is read only", call. = FALSE)
      }
    },

    #' @field varying alternative varying variables
    varying = function(value) {
      if (missing(value)) {
        private$.varying
      } else {
        stop("`verying` is read only", call. = FALSE)
      }
    }
  ),

  public = list(

    #' @field df data.frame containing the raw data
    df = data.frame(),

    #' @field alt_names strings containing the names of all alternatives.
    alt_names = list(),

    #' @field id string of the name of the column containing the decider ids.
    id = "",

    #' @field choice string of the name of the column containing the choices
    choice = "",

    #' @field ordered Boolean indicating if choices are ordered.
    ordered = FALSE,

    #' @description
    #' Create a new data object.
    #'
    #' @param df
    #' data.frame containing the raw data
    #' @param alt_names
    #' strings containing the names of all alternatives.
    #' @param id
    #' string of the name of the column containing the decider ids.
    #' @param choice
    #' string of the name of the column containing the choices
    #' @param ordered
    #' Boolean indicating if choices are ordered.
    #' @param varying
    #' strings of alternative varying variables
    #' @param dec_char
    #' strings of variables of decider characteristics
    initialize = function(
      df = data.frame(), alt_names = "", id = "", choice = "", ordered = FALSE,
      varying = "", dec_char = ""
    ) {
      stopifnot(
        is.data.frame(df), is.character(alt_names), is.character(id),
        is.character(choice), is.logical(ordered)
      )
      self$df <- df
      self$alt_names <- alt_names
      self$id <- id
      if (dim(df)[1] > 0) {
        ids <- df[, id]
        private$.N <- length(unique(ids))
        tt <- table(ids)
        private$.Tp <- as.vector(tt)
      } else {
        private$.N <- 0
        private$.Tp <- as.vector(0)
      }
      private$.varying <- varying
      private$.dec_char <- dec_char
      self$choice <- choice
      self$ordered <- ordered
    },

    #' @description
    #' Set df
    #' @param val data frame
    set_df = function(val) {
      if (is.data.frame(val)) {
        self$df <- val
        ids <- self$df[, self$id]
        private$.N <- length(unique(ids))
        tt <- table(ids)
        private$.Tp <- as.vector(tt)
      } else {
        cat("df must be a data.frame.")
      }
    },

    #' @description
    #' Set alt_names
    #' @param val alt_names
    set_alt_names = function(val) {
      if (is.list(val)) {
        self$alt_names <- val
      } else {
        cat("alt_names must be a string array.")
      }
    },

    #' @description
    #' Set id
    #' @param val id
    set_id = function(val) {
      if (is.character(val)) {
        self$id <- val
        ids <- self$df[, self$id]
        private$.N <- length(unique(ids))
        tt <- table(ids)
        private$.Tp <- as.vector(tt)
      } else {
        cat("id must be a string.")
      }
    },

    #' @description
    #' Set choice
    #' @param val choice
    set_choice = function(val) {
      if (is.character(val)) {
        self$choice <- val
        uc <- levels(self$df[, val])
        if (length(uc) > length(self$alt_names)) {
          self$alt_names <- uc
        }
      } else {
        cat("choice must be a string.")
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

    #' @description
    #' Set dec_char
    #' @param val dec_char
    set_dec_char = function(val) {
      if (is.character(val)) {
        private$.dec_char <- val
      } else {
        cat("dec_char must be a string.")
      }
    },

    #' @description
    #' Set varying
    #' @param val ordered
    set_varying = function(val) {
      if (is.character(val)) {
        private$.varying <- val
      } else {
        cat("varying must be a string.")
      }
    },

    #' @description
    #' print object
    print = function() {
      cat(sprintf(" N: %d, Tp: from %d to %d \n", private$.N, min(private$.Tp), max(private$.Tp)))
      cat(sprintf("Alternative varying variables: %s \n", paste(private$.varying, collapse = ", ")))
      cat(sprintf("Decider specific variables: %s \n", paste(private$.dec_char, collapse = ", ")))
      cat(sprintf("Alternatives: %s \n", paste(self$alt_names, collapse = ", ")))
      if (self$ordered) {
        cat(sprintf("Choice variables: %s (ordered)\n", self$choice))
      } else {
        cat(sprintf("Choice variables: %s \n", self$choice))
      }
      cat("==================================\n")
      print(summary(self$df))
    },

    #' @description
    #' present head of data frame object
    head = function() {
      head(self$df)
    },

    #' @description
    #' present tail of data frame object
    tail = function() {
      tail(self$df)
    },

    #' @description  plot the data frame: density plot for each chosen alternative.
    #' plot the data set
    plot = function() {
      choice <- self$df[, self$choice]
      if (is.factor(choice)) {
        alt <- levels(choice)
      } else {
        alt <- unique(choice)
      }
      # cycle over decider characteristics
      dec <- private$.dec_char
      if (dec != "") { # there are decider characteristics
        for (j in 1:length(dec)) {
          xch <- self$df[ind, dec[j]]
          ind <- which(choice == alt[1])
          if (var(xch[ind]) > 0) {
            plot(density(xch[ind]), col = 1, xlab = paste(dec[j]))
          } else {
            plot(xch[ind[1]], 0, col = 1, xlab = paste(dec[j]))
          }
          for (ja in 2:length(alt)) {
            ind <- which(choice == alt[ja])
            if (var(xch[ind]) > 0) {
              lines(density(xch[ind]), col = ja)
            }
          }
        }
      }
      vary <- private$.varying

      if (prod(vary == "") == 0) { # there are alternative varying regressors.
        bars <- rep(0, length(vary))
        tab <- list()
        for (jv in 1:length(vary)) {
          xv <- c()
          for (ja in 1:length(alt)) {
            xv <- c(xv, train_raw$df[, paste(vary[jv], "_", alt[ja], sep = "")])
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
              xv <- train_raw$df[, paste(vary[jv], "_", alt[1], sep = "")]
              if (var(xv[ind]) > 0) {
                plot(density(xv[ind]), col = 1, main = vary[jv], xlab = paste("Choice: ", alt[ja]))
              } else {
                plot(xv[ind[1]], 0, col = 1, main = vary[jv], xlab = paste("Choice: ", alt[ja]))
              }
              for (ja in 2:length(alt)) {
                xv <- train_raw$df[, paste(vary[jv], "_", alt[ja], sep = "")]
                if (var(xv[ind]) > 0) {
                  lines(density(xv[ind]), col = ja)
                }
              }
            } else {
              xv <- data.frame(x = train_raw$df[ind, paste(vary[jv], "_", alt[1], sep = "")], a = alt[1])
              for (ja in 2:length(alt)) {
                xv <- rbind(xv, data.frame(x = train_raw$df[, paste(vary[jv], "_", alt[ja], sep = "")], a = alt[ja]))
              }
              tt <- table(xv$x, xv$a)
              for (t in 1:length(alt)) {
                tt[, t] <- tt[, t] / sum(tt[, t])
              }
              barplot(t(tt), beside = TRUE, col = 1:length(alt))
            }
          }
        }
      }
      par(mfrow = c(1, 1))
    }
  )
)
