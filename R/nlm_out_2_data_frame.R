#' Transform nlm_output.txt file to data frame
#' @description Function that takes saved nlm output and extracts the information of the iterates and saves them in a data frame and if specified into a .csv file.
#' @param file_name name of the .txt file in which the nlm output is saved
#' @param save_in_file either a boolean, indicating if the data frame should be saved in a .csv-file or a character string indicating the name of the file in which the data frame should be saved.
#' @param return_dataframe boolean, indicating if the data frame should be returned.
#' @param parameter_names optional character vector, specifying the parameter names.
#' @return A data frame with iteration number, function value, parameter values and gradient values for each iteration of the nlm function.
#' @keywords internal

nlm_out_2_data_frame <- function(file_name, save_in_file = FALSE, return_dataframe = TRUE, parameter_names = NULL) {
  if (((save_in_file != FALSE) | (return_dataframe != FALSE)) & file.exists(file_name)) {
    raw_file <- utils::read.delim(file = file_name, header = FALSE)

    ### extract line by line
    data_frame <- c()
    i_iter <- 1
    while (i_iter <= nrow(raw_file)) {
      (text_line <- as.character(raw_file[i_iter, ]))

      if (substring(text_line, 1, 9) == "iteration") {
        which_step <- as.numeric(strsplit(x = text_line, split = " = ")[[1]][2])
      } else if (substring(text_line, 1, 9) == "Parameter") {
        ### initialise parameter vector
        Parameter_value <- c()

        ### go to next line in the code and extract line
        i_iter <- i_iter + 1
        text_line <- as.character(raw_file[i_iter, ])

        ### all lines till "Function Value" should belong to the parameter values
        while (substring(text_line, 1, 14) != "Function Value") {
          Parameter_value_add <- as.numeric(strsplit(x = text_line, split = " ")[[1]][-1])
          Parameter_value_add <- Parameter_value_add[!is.na(Parameter_value_add)]
          Parameter_value <- c(Parameter_value, Parameter_value_add)

          i_iter <- i_iter + 1
          text_line <- as.character(raw_file[i_iter, ])
        }
        i_iter <- i_iter - 1
      } else if (substring(text_line, 1, 14) == "Function Value") {
        i_iter <- i_iter + 1
        text_line <- as.character(raw_file[i_iter, ])
        Value <- as.numeric(strsplit(x = text_line, split = " ")[[1]][-1])
      } else if (substring(text_line, 1, 8) == "Gradient") {
        ### initialise gradient vector
        Gradient_value <- c()

        ### go to next line in the code and extract line
        i_iter <- i_iter + 1
        text_line <- as.character(raw_file[i_iter, ])

        ### all lines till "iteration = ..." should belong to the gradient
        while ((substring(text_line, 1, 9) != "iteration") & (i_iter <= nrow(raw_file))) {
          Gradient_value_add <- as.numeric(strsplit(x = text_line, split = " ")[[1]][-1])
          Gradient_value_add <- Gradient_value_add[!is.na(Gradient_value_add)]
          Gradient_value <- c(Gradient_value, Gradient_value_add)

          i_iter <- i_iter + 1
          text_line <- as.character(raw_file[i_iter, ])
        }

        ### check correctness of data
        length_parameter <- 0
        if (!is.null(data_frame)) {
          length_parameter <- (ncol(data_frame) - 2) / 2
        } else {
          length_parameter <- max(length(Parameter_value), length(Gradient_value))
        }

        if (length(Parameter_value) != length_parameter) {
          foo <- Parameter_value
          Parameter_value <- rep(NA, length_parameter)
          Parameter_value[1:min(length_parameter, length(foo))] <- foo[1:min(length_parameter, length(foo))]
          rm(foo)
        }
        if (length(Gradient_value) != length_parameter) {
          foo <- Gradient_value
          Gradient_value <- rep(NA, length_parameter)
          Gradient_value[1:min(length_parameter, length(foo))] <- foo[1:min(length_parameter, length(foo))]
          rm(foo)
        }


        data_frame_new <- c(which_step, Value, Parameter_value, Gradient_value)
        data_frame <- rbind(data_frame, data_frame_new)
        i_iter <- i_iter - 1
        which_step <- -1
        Value <- -1
        Parameter_value <- c()
        Gradient_value <- c()
      }

      (i_iter <- i_iter + 1)
    }


    ### set column and row names
    if (is.null(parameter_names)) {
      parameter_names <- 1:((ncol(data_frame) - 2) / 2)
    }
    colnames(data_frame) <- (c(
      "Step",
      "Value",
      paste("Parameter", parameter_names, sep = " "),
      paste("Gradient", parameter_names, sep = " ")
    ))
    rownames(data_frame) <- NULL
    data_frame <- as.data.frame(data_frame)


    if (is.character(save_in_file) | save_in_file == TRUE) {
      if (is.character(save_in_file)) {
        file_name_new <- save_in_file
      } else {
        file_name_new <- file_name
      }
      utils::write.csv(x = data_frame, file = paste(strsplit(file_name_new, split = ".txt")[[1]], ".csv", sep = ""), row.names = FALSE)
    }

    if (return_dataframe) {
      return(data_frame)
    }
  }
}
