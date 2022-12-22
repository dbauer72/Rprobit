#' Transform data for estimation
#' @description Function that transform 'data' by subtracting the regressor of the chosen alternative from the regressor matrix.
#' @param data \code{data}-object
#' @return transformed \code{data}-object \code{"data_tr"}
#' @keywords internal

substract_choice_regressor_from_data <- function(data) {
  N <- length(data)
  for (n in 1:N) {
    data_n <- data[[n]]
    T_n <- length(data_n$y)
    for (t in 1:T_n) {
      Xnt <- data_n$X[[t]]
      ynt <- as.numeric(data_n$y[t])
      Xntj <- Xnt[ynt, ]
      Xnt <- Xnt[-ynt, ] - matrix(1, dim(Xnt)[1] - 1, 1) %*% Xntj
      data_n$X[[t]] <- Xnt
    }
    data[[n]] <- data_n
  }
  return(data)
}
