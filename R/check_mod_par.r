#' Checks if parameterisation in the 'mod'-object is ok.
#' @description Function that checks if matrices in \code{mod}-object have correct dimensions.
#' @param mod \code{mod}-object
#' @return ok boolean indicating, if \code{mod}-object provides correct parameterisation
#' @keywords internal

check_mod_par <- function(mod = NULL) {

  ### check if dimensions fit ###
  alt <- mod$alt

  ok <- TRUE
  ### b: #rows in Hb and fB must be equal, #cols of Hb = lthb ###
  dimfb <- length(mod$fb)
  dimHb <- dim(mod$Hb)
  if (mod$lthb != dimHb[2]) {
    warning(paste0("lthb must match number of columns of Hb."))
    ok <- FALSE
  }
  if (dimfb != dimHb[1]) {
    warning(paste0("Length of fb must match number of rows of Hb."))
    ok <- FALSE
  }

  ### Omega: #rows in HO and fO must be equal and equal to lRE(lRE+1)/2, #cols of HO = lthO ###
  if (mod$lRE > 0) {
    dimfO <- length(mod$fO)
    dimHO <- dim(mod$HO)
    sizO <- mod$lRE * (mod$lRE + 1) / 2

    if (mod$lthO != dimHO[2]) {
      warning(paste0("lthO must match number of columns of HO."))
      ok <- FALSE
    }
    if (dimfO != dimHO[1]) {
      warning(paste0("Length of fO must match number of rows of HO."))
      ok <- FALSE
    }
    if (dimfO != sizO) {
      warning(paste0("Length of fO must match number of free entries for lRE random effects."))
      ok <- FALSE
    }
  }

  ### L: #rows in HL and fL must be equal, #cols of HL = lthL ###
  if (mod$ordered == FALSE) {
    dimfL <- length(mod$fL)
    dimHL <- dim(mod$HL)
    sizL <- alt * (alt + 1) / 2
    if (mod$lthL != dimHL[2]) {
      warning(paste0("lthL must match number of columns of HL."))
      ok <- FALSE
    }
    if (dimfL != dimHL[1]) {
      warning(paste0("Length of fL must match number of rows of HL."))
      ok <- FALSE
    }
    if (dimfL != sizL) {
      warning(paste0("Length of fL must match number of free entries for noise variance matrix."))
      ok <- FALSE
    }
  }

  if (ok == FALSE) {
    stop(paste0("Something is wrong with the model specification. See warnings for details. "))
  }

  return(ok)
}
