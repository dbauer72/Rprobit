#' Create mapping matrices
#' @description Function that creates mapping matrices Hb, fb, HO, fO, HL, fL if they do not exist in 'mod'.
#' These mappings define the structure of the Sigma matrix. Options are: \code{"fixed"} (\code{0.5*I}), \code{"ident"} (\code{sigma^2*I}), \code{"diag"}, \code{"full"}.
#' @param mod \code{mod}-object
#' @param corr_re Boolean, determining whether correlation for the random effects are estimated or not
#' @param error_struc Character, specifying the error structure
#' @return \code{mod}-object
#' @keywords internal

make_mappings <- function(mod, corr_re, error_struc) {

  ### identify missing values
  missing_values <- names(mod[is.na(mod)])

  ### setup Hb and fb
  if (any(c("Hb", "fb") %in% missing_values)) {
    if ("lthb" %in% missing_values) stop("Need to know the value of 'lthb' in 'mod'.")
    Hb <- diag(mod$lthb)
    fb <- rep(0, mod$lthb)
    mod$Hb <- Hb
    mod$fb <- fb
  }

  if (exists("ordered", where = mod) == FALSE) mod$ordered <- FALSE


  if (mod$ordered == TRUE) {
    nL <- mod$Tp[1]
  } else {
    if ("alt" %in% missing_values) stop("Need to know the value of 'alt' in 'mod'.")
    nL <- mod$alt
  }

  ### helper function: for a kxk matrix, it returns the position of the diagonal entries in the half vectorization
  fdiag <- function(k) {
    out <- rep(0, k)
    for (i in 1:k) out[i] <- (i - 1) * k + 1 - (i - 2) * (i - 1) / 2
    return(out)
  }

  ### setup HO, fO and lthO
  if (any(c("HO", "fO", "lthO") %in% missing_values)) {
    if ("lRE" %in% missing_values) stop("Need to know the value of 'lRE' in 'mod'.")
    if (mod$lRE == 0) {
      HO <- matrix(0, 0, 0)
      fO <- matrix(0, 0, 0)
    } else {
      HO <- diag((mod$lRE * mod$lRE - mod$lRE) / 2 + mod$lRE)
      if (corr_re == FALSE) {
        HO <- HO[, fdiag(mod$lRE), drop = FALSE]
      }
      fO <- rep(0, (mod$lRE * mod$lRE - mod$lRE) / 2 + mod$lRE)
    }
    lthO <- sum(HO, na.rm = TRUE)
    mod$HO <- HO
    mod$fO <- fO
    mod$lthO <- lthO
  }

  ### setup HL, fL and lthL
  if (any(c("HL", "fL", "lthL") %in% missing_values)) {
    if ("alt" %in% missing_values) stop("Need to know the value of 'alt' in 'mod'.")
    HL <- diag((nL * nL - nL) / 2 + nL)
    if (error_struc == "full") {
      HL <- HL[, -c(1:(nL + 1)), drop = FALSE]
      fL <- rep(0, (nL * nL - nL) / 2 + nL)
      fL[nL + 1] <- 1
      lthL <- sum(HL)
    }
    if (error_struc == "diag") {
      HL <- HL[, fdiag(nL), drop = FALSE]
      HL <- HL[, -1, drop = FALSE]
      fL <- rep(0, nrow(HL))
      # fL    <- rep(0,ncol(HL))
      fL[1] <- 1
      lthL <- sum(HL)
    }
    if (error_struc == "ident") {
      HL <- as.matrix(apply(HL[, fdiag(nL), drop = FALSE], 1, sum))
      # HL    <- HL[, -1, drop=FALSE]
      fL <- rep(0, nrow(HL))
      # fL    <- rep(0,ncol(HL))
      # fL[1] <- 1
      lthL <- 1
    }

    if (error_struc == "fixed") {
      HL <- HL[, fdiag(nL), drop = FALSE]
      fL <- HL %*% rep(sqrt(0.5), ncol(HL))
      HL <- matrix(0, (nL * nL - nL) / 2 + nL, 0)
      lthL <- 0
    }
    mod$HL <- HL
    mod$fL <- fL
    mod$lthL <- lthL
  }

  return(mod)
}
