#' Formatted printing of a matrix
#' @description Function that prints a formatted matrix.
#' @param est_sd matrix
#' @return inline written text
#' @keywords internal

print_est <- function(est_sd) {

  # check that `est_sd` has columns and rows
  if (any(dim(est_sd) == 0)) {
    cat("not available \n")
    return(NULL)
  }

  nr <- dim(est_sd)[1]
  nc <- dim(est_sd)[2]

  ### format the entries to have equal length
  est_sd_print <- est_sd
  for (r in 1:nr) {
    for (c in 1:(nc)) {
      est_sd_print[r, c] <- sprintf("%f", est_sd[r, c])
    }
  }
  for (c in 1:(nc)) {
    max_length <- max(nchar(est_sd_print[, c]))
    for (r in 1:nr) {
      while (nchar(est_sd_print[r, c]) < max_length) {
        est_sd_print[r, c] <- paste(" ", est_sd_print[r, c], sep = "")
      }
    }
  }

  ### print the matrix row by row
  for (r in 1:nr) {
    str <- sprintf("%d: ", r)
    for (c in 1:nc) {
      str <- c(str, sprintf(" %s ", est_sd_print[r, c]))
    }
    cat(c(str, "\n"))
  }
}
