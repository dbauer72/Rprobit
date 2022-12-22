#' Create identifier for decision maker
#' @description Function that combines multiple identifier in 'data_raw' into a single identifier.
#' @param data_raw \code{data_raw}-object
#' @param ids Vector, specifying the columns in 'data_raw' that identify a decision maker uniquely
#' @return vector, containing a unique identifier for each decision maker
#' @keywords internal

make_id <- function(data_raw, ids) {
  id <- 0
  for (i in seq_len(length(ids))) {
    id_add <- data_raw[ids[i]][, 1]
    if (!is.numeric(id_add)) {
      id_add <- as.numeric(as.factor(id_add))
    }
    id <- id * 10^(ceiling(log10(max(id_add)))) + id_add
  }
  return(id)
}
