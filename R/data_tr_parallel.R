#' Transform data for parallel framework
#' @description Function that transforms data for parallel framework.
#' @param Rprobit_obj \code{Rprobit_obj}-object
#' @param data_tr transformed \code{data}-object
#' @return A list containing the following elements:
#' \item{Rprobit_obj}{\code{Rprobit_obj}-object}
#' \item{data_tr_mat_list}{list of transformed data information}
#' @keywords internal

data_tr_parallel <- function(Rprobit_obj, data_tr) {

  ### change data to matrix form
  data_tr_X_mat <- matrix(0, 1, ncol(data_tr[[1]]$X[[1]]))
  data_tr_X_mat_indi <- matrix(0, Rprobit_obj$mod$N, 2)
  data_tr_X_mat_obs <- matrix(0, 1, 2)
  data_tr_y_mat <- matrix(0, 1, 1)
  data_tr_y_mat_indi <- matrix(0, Rprobit_obj$mod$N, 2)

  if (!is.null(Rprobit_obj$control$pairs_list)) {
    pairs_matrix <- matrix(0, 1, 3)
    pairs_matrix_indi <- matrix(0, Rprobit_obj$mod$N, 2)
  }

  for (i_N in 1:Rprobit_obj$mod$N) {
    data_tr_X_mat_indi[i_N, 1] <- nrow(data_tr_X_mat_obs)
    data_tr_y_mat_indi[i_N, 1] <- nrow(data_tr_y_mat)
    data_indi_N <- data_tr[[i_N]]

    if (!is.null(Rprobit_obj$control$pairs_list[[i_N]])) {
      pairs_matrix_iN <- Rprobit_obj$control$pairs_list[[i_N]]
      if (ncol(pairs_matrix_iN) == 3) {
        pairs_matrix_indi[i_N, 1] <- nrow(pairs_matrix)
        pairs_matrix <- rbind(pairs_matrix, pairs_matrix_iN)
        pairs_matrix_indi[i_N, 2] <- nrow(pairs_matrix) - 1
      } else {
        pairs_matrix_indi[i_N, ] <- 0
      }
    }

    for (i_obs in 1:length(data_indi_N$X)) {
      data_tr_X_mat_obs_N <- matrix(0, 1, 2)
      data_tr_X_mat_obs_N[1, 1] <- nrow(data_tr_X_mat)
      data_tr_X_mat <- rbind(data_tr_X_mat, data_indi_N$X[[i_obs]])
      data_tr_X_mat_obs_N[1, 2] <- nrow(data_tr_X_mat) - 1
      data_tr_X_mat_obs <- rbind(data_tr_X_mat_obs, data_tr_X_mat_obs_N)
    }
    data_tr_X_mat_indi[i_N, 2] <- nrow(data_tr_X_mat_obs) - 1
    data_tr_y_mat_indi[i_N, 1] <- nrow(data_tr_y_mat)
    data_tr_y_mat <- rbind(data_tr_y_mat, data_indi_N$y)
    data_tr_y_mat_indi[i_N, 2] <- nrow(data_tr_y_mat) - 1
  }
  data_tr_X_mat <- data_tr_X_mat[-1, ]
  data_tr_X_mat_obs <- data_tr_X_mat_obs[-1, ]
  data_tr_y_mat <- data_tr_y_mat[-1, ]

  if (!is.null(Rprobit_obj$control$pairs_list)) {
    pairs_matrix <- pairs_matrix[-1, ]
  }

  data_tr_mat_list <- list(
    "data_tr_X_mat" = data_tr_X_mat,
    "data_tr_X_mat_indi" = data_tr_X_mat_indi,
    "data_tr_X_mat_obs" = data_tr_X_mat_obs,
    "data_tr_y_mat" = data_tr_y_mat,
    "data_tr_y_mat_indi" = data_tr_y_mat_indi
  )

  if (!is.null(Rprobit_obj$control$pairs_list)) {
    Rprobit_obj$control$pairs_list <- list(
      "pairs_matrix" = pairs_matrix,
      "pairs_matrix_indi" = pairs_matrix_indi
    )
  }

  return(list("Rprobit_obj = Rprobit_obj", "data_tr_mat_list" = data_tr_mat_list))
}
