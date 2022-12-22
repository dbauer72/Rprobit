#' Calculate Confusion Matrix of Model Predictions
#' 
#' @description 
#' This function computes the confusion matrix comparing the predicted choices 
#' with the actual choices.
#' 
#' @param Rprobit_obj 
#' An \code{Rprobit_cl} object
#' @param data_new 
#' Optionally new data for prediction. Can be
#' * a \code{\link{data_cl}} object
#' * or a \code{\link{data_raw_cl}} object.
#' By default, \code{data_new = NULL}.
#' 
#' @return 
#' A \code{list} with elements:
#' * \code{pred_choices}, a \code{numeric} vector of predicted choices,
#' * \code{choices}, a \code{numeric} vector of actual choices,
#' * \code{TT}, a \code{table}, the confusion matrix.
#' 
#' @export

confmat_Rprobit <- function(Rprobit_obj, data_new = NULL){

  ### make sure the Rprobit_obj object is a "Rprobit_cl" object
  if(!inherits(Rprobit_obj, "Rprobit_cl")){
    stop("'Rprobit_obj' is not a Rprobit object", call. = FALSE)
  }
  
  ### clone object to avoid side effects.
  Rprobit_o <- Rprobit_obj$clone(deep=TRUE)
  
  ### check if list version of data is available:
  data_from_data_raw_model <- function(Rprobit_o){
    read_formula_out  <- read_formula(Rprobit_o$form)
    ASC               <- read_formula_out$ASC
    allvars           <- read_formula_out$allvars
    choice            <- all.vars(Rprobit_obj$form)[1]
    norm_alt          <- Rprobit_obj$info$setup_input$norm_alt
    alt               <- Rprobit_obj$mod$alt
    data_raw_to_data(
      data_raw = Rprobit_o$data_raw, allvars = allvars, choice = choice, 
      re = Rprobit_o$re, norm_alt = norm_alt, alt = alt, ASC = ASC
    )
  }

  ### add new data
  if (!(is.null(data_new))) {
    if (inherits(data_new, "data_cl")) {
      Rprobit_o$data <- data_new
    }
    if (inherits(data_new, "data_raw_cl")){
      Rprobit_o$data_raw  = data_new
      Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
    }
  } else if(is.null(Rprobit_o[["data"]]) & !is.null(Rprobit_o[["data_raw"]])){
    Rprobit_o$data <- data_from_data_raw_model(Rprobit_o)
  }

  ### predict
  prediction_probs <- matrix(NA, 0, Rprobit_o$mod$alt)
  choices <- numeric(0)
  for (n in 1:Rprobit_o$data$N){

    ### calculate predictions for each choice occasion
    data_n <- Rprobit_o$data$data[[n]]

    if (Rprobit_obj$mod$ordered == FALSE){
      pred_n <- pred_probit_approx(Rprobit_o$theta,data_n,Rprobit_o$mod,Rprobit_o$control$approx_method)
      prediction_probs <- rbind(prediction_probs,pred_n)
    } else {
      pred_n <- pred_probit_ordered_approx(Rprobit_o$theta,data_n$X,data_n$y,Rprobit_o$mod)
      prediction_probs <- rbind(prediction_probs,pred_n)
    }

    ### store true choices
    choices <- c(choices, data_n$y)
  }

  ### calculate the predicted choices
  pred_choices <- numeric(0)
  for (j in 1:dim(prediction_probs)[1]) {
    pred_choices <- c(pred_choices, which(prediction_probs[j, ] == max(prediction_probs[j, ])))
  }
  TT <- table(pred_choices,choices)
  colnames(TT) <- Rprobit_obj$data_raw$alt_names[as.numeric(colnames(TT))]
  rownames(TT) <- Rprobit_obj$data_raw$alt_names[as.numeric(rownames(TT))]
  
  ### return true choice percentages and averages of prediction
  list(pred_choices = pred_choices, choices = choices, conf_mat = TT)
}
