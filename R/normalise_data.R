#' Normalise data
#' @description Function that normalises data_raw regressor variables
#' @param data_raw \code{data_raw_cl}-object
#' @param normalisation \code{normalisation}-object. Either a list with normalisation specifications or a boolean.
#' @param allvars \code{allvars}-object
#' @param ASC boolean to check if ASCs are included in the model
#' @return \code{data_raw_cl}-object with normalisation specifications as attribute.
#' @keywords internal

normalise_data <- function(data_raw, normalisation = FALSE, allvars = NULL, ASC = TRUE) {

  ### check if data has already been normalised
  check_normalised <- TRUE
  if (!is.null(attr(x = data_raw, which = "normalised"))) {
    if (attr(x = data_raw, which = "normalised")) {
      check_normalised <- FALSE
    }
  }


  ### normalise data
  if ((normalisation == TRUE) & check_normalised) {

    ### warning if no ASCs are present but notmalisation is used
    if (!ASC & check_normalised) {
      message("Warning: Normalising data without including ASCs changes the model by base shifting utilities.")
    }

    if (is.list(normalisation)) {
      ### if a list of normalisations is given

      for (i_var in names(normalisation_list)[names(normalisation_list) %in% colnames(data_raw$df)]) {
        print(i_var)
        i_col <- which(colnames(data_raw$df) == i_var)

        a_norm <- normalisation[[i_var]]$a_norm
        b_norm <- normalisation[[i_var]]$b_norm
        if (is.na(a_norm)) a_norm <- 0
        if (b_norm == 0 | is.na(b_norm)) b_norm <- 1
        normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
        data_raw$df[, i_col] <- (data_raw$df[, i_col] - a_norm) / b_norm
      }

      for (i_var in colnames(data_raw$df)[!(colnames(data_raw$df) %in% names(normalisation_list))]) {
        normalisation_list[[i_var]] <- list(a_norm = 0, b_norm = 1)
      }

      attr(x = data_raw, which = "normalisation") <- normalisation_list
      attr(x = data_raw, which = "normalised") <- TRUE
    } else if (is.character(normalisation)) {

      ### for special normailsations (not implemented)
    } else if (normalisation) {
      ### if normalisation==TRUE

      normalisation_list <- list()
      for (i_var in allvars[[1]]) {
        if (i_var != 0) {
          norm_cols <- grep(paste0(i_var, "_"), colnames(data_raw$df))
          if (any(!(unique(data_raw$df[, norm_cols]) %in% c(0, 1)))) {

            ### use overall normalisation for variables of type 1 (since only one regressor is estimated)
            a_norm <- mean(as.matrix(data_raw$df[, norm_cols]))
            b_norm <- stats::sd(as.matrix(data_raw$df[, norm_cols]))
            if (is.na(a_norm)) a_norm <- 0
            if (b_norm == 0 | is.na(b_norm)) b_norm <- 1

            for (i_col in norm_cols) {
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
              data_raw$df[, i_col] <- (data_raw$df[, i_col] - a_norm) / b_norm
            }
          } else {
            for (i_col in norm_cols) {
              a_norm <- 0
              b_norm <- 1
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
            }
          }
        }
      }

      for (i_var in allvars[[2]]) {
        if (i_var != 0) {
          norm_cols <- grep(i_var, colnames(data_raw$df))
          for (i_col in norm_cols) {
            if (any(!(unique(data_raw$df[, i_col]) %in% c(0, 1)))) {
              a_norm <- mean(data_raw$df[, i_col])
              b_norm <- stats::sd(data_raw$df[, i_col])
              if (is.na(a_norm)) a_norm <- 0
              if (b_norm == 0 | is.na(b_norm)) b_norm <- 1
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
              data_raw$df[, i_col] <- (data_raw$df[, i_col] - a_norm) / b_norm
            } else {
              a_norm <- 0
              b_norm <- 1
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
            }
          }
        }
      }

      for (i_var in allvars[[3]]) {
        if (i_var != 0) {
          norm_cols <- grep(paste0(i_var, "_"), colnames(data_raw$df))
          for (i_col in norm_cols) {
            if (any(!(unique(data_raw$df[, i_col]) %in% c(0, 1)))) {
              a_norm <- mean(data_raw$df[, i_col])
              b_norm <- stats::sd(data_raw$df[, i_col])
              if (is.na(a_norm)) a_norm <- 0
              if (b_norm == 0 | is.na(b_norm)) b_norm <- 1
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
              data_raw[, i_col] <- (data_raw$df[, i_col] - a_norm) / b_norm
            } else {
              a_norm <- 0
              b_norm <- 1
              normalisation_list[[colnames(data_raw$df)[i_col]]] <- list(a_norm = a_norm, b_norm = b_norm)
            }
          }
        }
      }

      attr(x = data_raw, which = "normalisation") <- normalisation_list
      attr(x = data_raw, which = "normalised") <- TRUE
    }
  } else {
    attr(x = data_raw, which = "normalisation") <- FALSE
  }


  return(data_raw)
}
