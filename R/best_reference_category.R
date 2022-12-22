#' Find Best Reference Categories for Categorical Variables
#'
#' @description
#' This function finds best reference categories for categorical variables,
#' optimally a combination of categories that has optimal coverage over all
#' choice alternatives. To be used before creating dummy variables!
#'
#' @param data_raw
#' A \code{\link{data_raw_cl}} object.
#' @param categorical_variables
#' A \code{character} vector, containing the names of the categorical variables.
#' @param choice_variable
#' A \code{character}, the name of the choice variable.
#' @param batches
#' A \code{list} of boolean vectors to specify subsets of the raw data.
#' By default, \code{batches = NULL}.
#'
#' @return
#' A \code{matrix} containing information on the number of choices made for each
#' combination of categories and each choice alternative, ordered by recommendation.
#'
#' @keywords internal

best_reference_category <- function(
    data_raw, categorical_variables, choice_variable, batches = NULL
  ) {


  ### choice alternatives
  choice_alternatives <- sort(unique(data_raw[, choice_variable]))

  ### categorical variables
  category_list <- list()
  for (i1 in 1:length(categorical_variables)) {
    category_list[[i1]] <- sort(unique(data_raw[, categorical_variables[i1]]))
  }

  ### build matrix
  choice_overview <- expand.grid(category_list)
  choice_overview <- cbind(choice_overview, matrix(0, nrow(choice_overview), length(choice_alternatives) + 2))
  if (is.numeric(categorical_variables)) {
    colnames(choice_overview) <- c(colnames(data_raw)[categorical_variables], as.character(choice_alternatives), "total", "n_missing")
  } else if (is.character(categorical_variables)) {
    colnames(choice_overview) <- c(categorical_variables, as.character(choice_alternatives), "total", "n_missing")
  }


  ### check choices
  for (i1 in 1:nrow(choice_overview)) {

    ### select according subset
    subset_row <- rep(TRUE, nrow(data_raw))
    for (i2 in 1:length(categorical_variables)) {
      subset_row <- subset_row & (data_raw[, categorical_variables[i2]] == choice_overview[i1, i2])
    }

    ### count number of choices for each alternative
    for (i2 in 1:length(choice_alternatives)) {
      choice_overview[i1, length(categorical_variables) + i2] <- sum(data_raw[subset_row, choice_variable] == choice_alternatives[i2])
    }
    choice_overview[i1, ncol(choice_overview) - 1] <- sum(subset_row)
    choice_overview[i1, ncol(choice_overview)] <- sum(choice_overview[i1, (length(categorical_variables) + 1):(length(categorical_variables) + length(choice_alternatives))] == 0)

    if (is.null(batches)) {
      message(sprintf("Completed: %.0f%%", (i1 / nrow(choice_overview) * 100)), "\r", appendLF = FALSE)
    } else {
      message(sprintf("Completed: %.0f%%", (i1 / nrow(choice_overview) * 100 / (1 + length(batches)))), "\r", appendLF = FALSE)
    }
  }


  ### some more metrics
  average_cat_freq <- numeric(nrow(choice_overview))
  for (i in 1:length(choice_alternatives)) {
    average_cat_freq <- average_cat_freq + choice_overview[, length(categorical_variables) + i] / sum(choice_overview[, length(categorical_variables) + i]) / length(choice_alternatives)
  }
  average_log_cat_freq <- numeric(nrow(choice_overview))
  for (i in 1:length(choice_alternatives)) {
    average_log_cat_freq <- average_log_cat_freq + log(choice_overview[, length(categorical_variables) + i] / sum(choice_overview[, length(categorical_variables) + i])) / length(choice_alternatives)
  }
  choice_overview <- cbind(choice_overview, average_cat_freq, average_log_cat_freq)


  ### sort new
  choice_overview <- choice_overview[order(choice_overview$total, decreasing = TRUE), ]
  choice_overview <- choice_overview[order(choice_overview$n_missing, decreasing = FALSE), ]
  choice_overview <- choice_overview[order(choice_overview$average_log_cat_freq, decreasing = TRUE), ]
  rownames(choice_overview) <- NULL


  ### if you have batches
  n_missing_in_batches <- numeric(nrow(choice_overview))
  average_cat_freq_batches <- numeric(nrow(choice_overview))
  average_log_cat_freq_batches <- numeric(nrow(choice_overview))
  if (!is.null(batches)) {
    for (i_batch in 1:length(batches)) {
      choice_overview_batch <- choice_overview

      for (i1 in 1:nrow(choice_overview_batch)) {
        subset_row <- batches[[i_batch]]
        for (i2 in 1:length(categorical_variables)) {
          subset_row <- subset_row & (data_raw[, categorical_variables[i2]] == choice_overview_batch[i1, i2])
        }

        for (i2 in 1:length(choice_alternatives)) {
          choice_overview_batch[i1, length(categorical_variables) + i2] <- sum(data_raw[subset_row, choice_variable] == choice_alternatives[i2])
        }

        choice_overview_batch$total[i1] <- sum(subset_row)
        choice_overview_batch$n_missing[i1] <- sum(choice_overview_batch[i1, (length(categorical_variables) + 1):(length(categorical_variables) + length(choice_alternatives))] == 0)

        message(sprintf("Completed: %.0f%%", (100 * (i_batch + i1 / nrow(choice_overview)) / (1 + length(batches)))), "\r", appendLF = FALSE)
      }

      n_missing_in_batches <- n_missing_in_batches + choice_overview_batch$n_missing

      ### some more metrics
      for (i in 1:length(choice_alternatives)) {
        average_cat_freq_batches <- average_cat_freq_batches + choice_overview_batch[, length(categorical_variables) + i] / sum(choice_overview_batch[, length(categorical_variables) + i]) / length(choice_alternatives)
      }
      for (i in 1:length(choice_alternatives)) {
        average_log_cat_freq_batches <- average_log_cat_freq_batches + log(choice_overview_batch[, length(categorical_variables) + i] / sum(choice_overview_batch[, length(categorical_variables) + i])) / length(choice_alternatives)
      }
    }
    choice_overview <- cbind(choice_overview, n_missing_in_batches, "average_cat_freq_batches" = average_cat_freq_batches / length(batches), "average_log_cat_freq_batches" = average_log_cat_freq_batches / length(batches))

    ### reorder
    choice_overview <- choice_overview[order(choice_overview$n_missing_in_batches, decreasing = FALSE), ]
    choice_overview <- choice_overview[order(choice_overview$average_log_cat_freq, decreasing = TRUE), ]
    choice_overview <- choice_overview[order(choice_overview$average_log_cat_freq_batches, decreasing = TRUE), ]
    rownames(choice_overview) <- NULL
  }


  return(choice_overview)
}
