#' Convert Long Data
#'
#' @description
#' This function converts data in long format to an object of class
#' \code{data_raw_cl}.
#'
#' @param df
#' A \code{data.frame} in long format.
#' By default, \code{df = NULL}.
#' @param id_dec
#' A \code{character}, the column name in \code{df} containing the ids uniquely
#' identifying each decider.
#' By default, \code{id_dec = "id_dec"}.
#' @param id_choice
#' A \code{character}, the column name in \code{df} identifying the choice situations.
#' By default, \code{id_choice = "id_choice"}.
#' @param alternative
#' A \code{character}, the column name in \code{df} identifying the various
#' alternatives in the long format.
#' By default, \code{alternative = "mode"}.
#' @param choice
#' A \code{character}, the column name in \code{df} column containing the choices.
#' By default, \code{choice = "choice"}.
#'
#' @return
#' A \code{\link{data_raw_cl}} object.
#'
#' @export

convert_data_long_data_raw <- function(
    df = NULL, id_dec = "id_dec", id_choice = "id_choice", alternative = "mode",
    choice = "choice")
  {

  if (is.null(df)) {
    stop("Data must be provided.", call. = FALSE)
  }

  ids <- df[, id_dec]
  ids_sit <- df[, id_choice]
  alt_names <- unique(df[, alternative]) # alternatives
  if (is.factor(alt_names)) {
    alt_names <- levels(alt_names)
  }

  individuals <- unique(ids)
  id_joint <- c()
  for (j in 1:length(individuals)) {
    ind_ind <- which(df[, id_dec] == individuals[j])
    situations <- unique(df[ind_ind, id_choice])
    for (jj in 1:length(situations)) {
      id_joint <- rbind(id_joint, c(individuals[j], situations[jj]))
    }
  }
  df_raw <- data.frame(id = id_joint[, 1], id_sit = id_joint[, 2])
  # cycle over variables
  vars <- colnames(df)
  out <- c(id_choice, id_dec, choice, alternative)
  varying <- c() # names of alternative varying variables
  dec_char <- c() # names of decider characteristic variables.

  alt <- length(alt_names)
  for (j in 1:length(vars)) {
    if ((vars[j] %in% out) == FALSE) {
      vals <- matrix(df[, vars[j]], nrow = alt)
      # check for variation
      iota <- rep(1, alt)
      vals_dm <- (diag(alt) - iota %*% t(iota) / alt) %*% vals
      if (sum(abs(vals_dm)) > 0) { # then the variable varies over alternatives
        varying <- c(varying, vars[j])
        for (a in 1:alt) {
          df_raw[paste0(vars[j], "_", alt_names[a])] <- vals[a, ]
        }
      } else { # then it is a decider characteristic.
        dec_char <- c(dec_char, vars[j])
        df_raw[vars[j]] <- vals[1, ]
      }
    }
  }

  # last bit: include decision.
  df_raw["choice"] <- 0
  for (j in 1:dim(df_raw)[1]) {
    ind_ind <- which(df[, id_choice] == df_raw[j, "id_sit"])
    ch <- which(df[ind_ind, choice] == "yes")
    df_raw[j, "choice"] <- alt_names[as.numeric(df[ind_ind[ch], alternative])]
  }
  df_raw[, "choice"] <- as.factor(df_raw[, "choice"])

  # now generate object
  data_raw_cl$new(
    df = df_raw,
    alt_names = alt_names,
    id = "id",
    choice = "choice",
    ordered = FALSE,
    varying = varying,
    dec_char = dec_char
  )
}
