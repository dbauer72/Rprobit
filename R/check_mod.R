#' Checks if input is a valid 'mod'-object
#' @description Function that checks if 'mod' is a valid \code{mod}-object.
#' @param mod \code{mod}-object
#' @return (valid) \code{mod}-object
#' @keywords internal

check_mod <- function(mod = NULL) {
  ### ingredients of 'mod'
  ingredients <- c("alt", "lthb", "lthO", "lthL", "lRE", "Hb", "fb", "HO", "fO", "HL", "fL", "ordered")

  ### create blueprint of 'mod'
  if (is.null(mod)) {
    mod <- as.list(rep(NA, length(ingredients)))
    names(mod) <- ingredients
  } else {
    ### add "ordered" if not there.
    if (exists("ordered", where = mod) == FALSE) mod$ordered <- FALSE
    ### add missing elements of 'mod'
    if (length(setdiff(ingredients, names(mod))) > 0) {
      warning(paste0("The following elements are missing in 'mod': ", paste(setdiff(ingredients, names(mod)), collapse = " "), "."))
      mod[setdiff(ingredients, names(mod))] <- NA
    }

    ### remove unsupported elements of 'mod'
    if (length(setdiff(names(mod), ingredients)) > 0) {
      warning(paste0("The following elements in 'mod' are not supported and will be removed: ", paste(setdiff(names(mod), ingredients), collapse = " "), "."))
      mod[setdiff(names(mod), ingredients)] <- NULL
    }
  }
  if (exists("ordered", where = mod) == TRUE) mod$ordered <- FALSE
  return(mod)
}
