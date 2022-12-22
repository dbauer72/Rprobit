#' Check availability of alternatives
#' @description Function that checks the availability of alternatives for each decision maker and choice occasion.
#' @param data_raw \code{data_raw}-object
#' @param avail List, contains information on the availability for each alternative and each choice occasion
#' @param Rprobit_obj \code{Rprobit_obj}-object
#' @return Rprobit_obj
#' @keywords internal

check_avail <- function(data_raw, avail, Rprobit_obj) {

  ### preperations for NA if not available
  if (!is.null(avail)) {

    ### check if this is a list
    if (is.list(avail)) {

      ### check if there is a list entry for every alternative
      if (length(avail) == Rprobit_obj$mod$alt) {
        for (i_alt in 1:Rprobit_obj$mod$alt) {

          ### get total number of choice occasions
          data_length <- sum(data_raw$Tp)
          if (length(data_raw$Tp) == 1) {
            data_length <- data_length * data_raw$N
          }

          ### check length, if not correct, assume alternative is available
          if (length(avail[[i_alt]]) != data_length) {
            avail[[i_alt]] <- rep(1, data_length)
          }
          ### if entry is 0, set to not available
          avail[[i_alt]][avail[[i_alt]] == 0] <- NA

          ## if entry is available, set to 1
          avail[[i_alt]][!is.na(avail[[i_alt]])] <- 1
        }

        ### go through all decision makers
        if (!is.null(data_raw)) {
          id <- data_raw$id
          id_data <- data_raw$df[, id]
          id_macml <- unique(id_data)
        } else {
          id_macml <- 1:length(Rprobit_obj$data)
          id_data <- c()
          for (i_id in 1:length(Rprobit_obj$data)) {
            id_data <- c(id_data, rep(i_id, length(Rprobit_obj$data[[i_id]]$X)))
          }
        }

        for (i_N in 1:data_raw$N) {
          Tp_N <- data_raw$Tp[min(length(data_raw$Tp), i_N)]

          ### go through all time points
          for (i_tp in 1:Tp_N) {

            ### go through all alternatives
            for (i_alt in 1:Rprobit_obj$mod$alt) {
              ### if it is not the reference alternative
              if (i_alt != Rprobit_obj$info$setup_input$norm_alt) {
                ### set alternative covariates to NA if it is not available
                avail_i_N <- (avail[[i_alt]])[id_data == id_macml[i_N]][i_tp]
                Rprobit_obj$data[[i_N]]$X[[i_tp]][i_alt, ] <- Rprobit_obj$data[[i_N]]$X[[i_tp]][i_alt, ] * avail_i_N
              }
            }
          }
        }
      }
    }
  }

  return(Rprobit_obj)
}
