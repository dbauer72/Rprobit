#' Calculate Choice Predictions
#' 
#' @description 
#' This function calculates and returns true choice percentages and averages of 
#' choice prediction.
#' 
#' @param Rprobit_obj 
#' An \code{\link{Rprobit_cl}} object.
#' @param data_new 
#' Optionally new data for prediction. Can be
#' * a \code{\link{data_cl}} object
#' * or a \code{\link{data_raw_cl}} object.
#' By default, \code{data_new = NULL}.
#' @param all_pred 
#' A \code{logical}, indicating whether mean predicted probabilities or the 
#' matrix of all predicted probabilities should be returned.
#' 
#' @return 
#' A \code{matrix}, true choice percentages in first and averages of predictions 
#' in second row
#' 
#' @export

predict_Rprobit <- function(Rprobit_obj, data_new = NULL, all_pred = FALSE) {

  ### make sure the Rprobit_obj object is a "Rprobit" object
  if (!inherits(Rprobit_obj, "Rprobit_cl")) {
    stop("'Rprobit_obj' is not a Rprobit_cl object")
  }

  # clone object to avoid side effects.
  Rprobit_o <- Rprobit_obj$clone()

  ### check if list version of data is available:
  data_from_data_raw_model <- function(Rprobit_o) {
    read_formula_out <- read_formula(Rprobit_o$form)
    ASC <- read_formula_out$ASC
    allvars <- read_formula_out$allvars
    choice <- all.vars(Rprobit_obj$form)[1]
    norm_alt <- Rprobit_obj$info$setup_input$norm_alt
    alt <- Rprobit_obj$mod$alt

    data_raw_to_data_out <- data_raw_to_data(data_raw = Rprobit_o$data_raw, allvars = allvars, choice = choice, re = Rprobit_o$re, norm_alt = norm_alt, alt = alt, ASC = ASC)
    # data_raw_to_data_out  <- data_raw_to_data(data_raw = Rprobit_obj$data_raw, allvars = allvars, choice = choice, re = Rprobit_obj$re, norm_alt = norm_alt, alt = alt, ASC = ASC)

    return(data_raw_to_data_out)
  }


  ### add new data
  if (!(is.null(data_new))) {
    #    if (class(data_new) == "data_cl"){
    if (inherits(data_new, "data_cl")) {
      Rprobit_o$data <- data_new
      # Rprobit_obj$mod$Tp <- data_new$Tp
    }
    #    if(!(is.matrix(data_new) | is.data.frame(data_new)) & is.list(data_new) ){
    #      ### if it is already transformed data, directly add it to the model object
    #      Rprobit_obj$data  <- data_new
    #      Rprobit_obj$mod$N <- length(data_new)
    #    } else if(is.matrix(data_new) | is.data.frame(data_new)) {
    ### if it is raw data, transform it to the required format
    #    if (class(data_new) == "data_raw_cl"){


    if (inherits(data_new, "data_raw_cl")) {
      Rprobit_o$data_raw <- data_new
      Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
    }
  } else if (is.null(Rprobit_o[["data"]]) & !is.null(Rprobit_o[["data_raw"]])) {
    if (Rprobit_o$mod$ordered) {
      data_raw_obj <- Rprobit_o$data_raw$clone()
      data <- data_cl$new(ordered = TRUE, vars = data_raw_obj$dec_char)
      ids <- data_raw_obj$df[, data_raw_obj$id]
      unique_ids <- unique(ids)
      data_df <- list()
      for (j in 1:length(unique_ids)) {
        ind <- which(ids == unique_ids[j])
        data_df[[j]] <- list(X = data.matrix(data_raw_obj$df[ind, data_raw_obj$dec_char]), y = data.matrix(data_raw_obj$df[ind, data_raw_obj$choice]))
      }
      # vars = allvars[[1]]
      # alt_names = NULL
      data$set_data(data_df)
      Rprobit_o$data <- data
    } else {
      Rprobit_o$data <- data_from_data_raw_model(Rprobit_obj)
    }


    # Rprobit_obj$mod$Tp = Rprobit_obj$data$Tp
  }

  ### predict
  predictions <- matrix(NA, 0, Rprobit_o$mod$alt)
  choices <- numeric(0)
  for (n in 1:length(Rprobit_o$data$Tp)) {

    ### calculate predictions for each choice occasion
    data_n <- Rprobit_o$data$data[[n]]

    if (Rprobit_o$mod$ordered == FALSE) {
      pred_n <- pred_probit_approx(Rprobit_o$theta, data_n, Rprobit_o$mod, Rprobit_o$control$approx_method)
      predictions <- rbind(predictions, pred_n)
    } else {
      ### TODO: This function is not defined in this branch. To Dietmar: is it somewhere in one of your branches?
      pred_n <- pred_probit_ordered_approx(Rprobit_o$theta, data_n$X, data_n$y, Rprobit_o$mod)

      predictions <- rbind(predictions, pred_n)
    }

    ### store true choices
    choices <- c(choices, data_n$y)
  }

  ### return true choice percentages and averages of prediction
  if (all_pred) {
    out <- predictions
    colnames(out) <- Rprobit_o$data_raw$alt_names
    rownames(out) <- NULL
    out <- as.data.frame(out)
    out$choice <- Rprobit_o$data_raw$alt_names[choices]
    out$choice_prob <- as.numeric(out[cbind(1:length(choices), choices)])
  } else {
    out <- rbind(table(choices) / sum(table(choices)), colMeans(predictions))
    colnames(out) <- paste("alternative", colnames(out))
    colnames(out) <- Rprobit_o$data_raw$alt_names
    rownames(out) <- c("true choice percentages", "averages of predictions")
  }
  return(out)
}
