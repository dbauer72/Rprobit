#' Transform 'data_raw' to 'data'
#' @description Function that takes 'data_raw' and transforms it into 'data'. Operations include expansion of variables into regressors using suitable differencing and omissions to prevent dummy variable trap. Data_raw is assumed to be in wide format, with alternative specific variables indicated using 'var_cat' notation.
#' @param data_raw \code{data_raw_cl}-object
#' @param allvars \code{allvars}-object
#' @param choice Character name of choice variable
#' @param re Names of variables to be considered as random effects
#' @param norm_alt Norm alternative: either an integer or a string
#' @param alt Number of choice alternatives
#' @param ASC Boolean determining whether ASCs are estimated or not
#' @return \code{data_cl} object
#' @keywords internal

data_raw_to_data <- function(data_raw, allvars, choice, re, norm_alt, alt, ASC) {

  ### auxiliary function: create k x k identity matrix, cutting out the j-th column
  diagn <- function(k, j) if (j != 0) diag(k)[, -j] else diag(k)

  ### auxiliary function: create k x k identity matrix, replacing the j-th column by -1's and the j-th row by 0's
  diagdiff <- function(k, j) {
    if (j != 0) {
      out <- diag(k)
      out[, j] <- (-1)
      out[j, ] <- 0
      return(out)
    } else {
      diag(k)
    }
  }

  ### transform choice variable to numeric
  #choice_levels <- levels(as.factor(data_raw$alt_names))
  choice_levels <- data_raw$alt_names
  #choice_levels <- reorder(choice_levels, sort(as.numeric(choice_levels)))
  if(length(choice_levels)==0) choice_levels <- as.factor(1:alt)
  if(length(choice_levels)!=alt) choice_levels <- as.factor(1:alt)
  #data_raw[,choice] <- as.numeric(as.factor(data_raw[,choice]))
  
  ### extract meta data
  # id_macml <- unique(data_raw$id_macml)
  id <- data_raw$id
  id_macml <- data_raw$df[, id]
  ids <- unique(id_macml)
  # N <- length(id_macml)
  N <- data_raw$N
  # Tp <- numeric(N)
  Tp <- data_raw$Tp
  # if(is.numeric(data_raw[,choice])){
  #  alt_names <- levels(1:alt
  # } else {
  #  alt_names <- levels(data_raw[,choice])
  # }
  alt_names <- choice_levels

  # convert norm_alt to numeric, if necessary
  if (!is.numeric(norm_alt)) {
    norm_alt <- which(alt_names == norm_alt)
  }

  if (norm_alt != 0) alt_names <- alt_names[-norm_alt]
  data <- list()

  data_raw_n <- data_raw$df[data_raw$df[, id] == ids[1], ]
  Xnt <- data_raw_n[1, ]

  ord <- list()
  ord[["ASC"]] <- c()
  for (j in alt_names) {
    ord[["ASC"]] <- c(ord[["ASC"]], grep(paste0("ASC", "_", j), colnames(Xnt)))
  }

  ### find ordering of alternative specific variables
  for (i_var in allvars[[1]]) {
    if (i_var != 0) {
      ord[[i_var]] <- c()
      for (j in choice_levels) {
        ord[[i_var]] <- c(ord[[i_var]], grep(paste0(i_var, "_", j), colnames(Xnt)))
      }
    }
  }

  for (i_var in allvars[[3]]) {
    if (i_var != 0) {
      ord[[i_var]] <- c()
      for (j in choice_levels) {
        ord[[i_var]] <- c(ord[[i_var]], grep(paste0(i_var, "_", j), colnames(Xnt)))
      }
    }
  }
  ### loop over decision makers
  for (n in 1:N) {
    data_raw_n <- data_raw$df[id_macml == ids[n], ]
    Tp[n] <- dim(data_raw_n)[1]
    yn <- data_raw_n[, choice]
    Xn <- list()

    ### loop over choice occasions
    for (t in 1:Tp[n]) {
      Xnt <- data_raw_n[t, ]
      Xnt_temp <- data.frame(remove_me = 1:alt)

      ### add all the alternative specific variables
      for (i_var in allvars[[1]]) {
        if (i_var != 0) {
          old_col <- colnames(Xnt_temp)
          ### locate the relevant columns and compute difference w.r.t norm_alt, put covariates with random effects at the end
          if (i_var %in% re) {
            # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
            Xnt_temp <- cbind(diagdiff(alt, norm_alt) %*% t(Xnt[, ord[[i_var]], drop = FALSE]), Xnt_temp)
            # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var,"_"), colnames(Xnt))])), Xnt_temp)
            colnames(Xnt_temp) <- c(i_var, old_col)
          } else {
            # Xnt_temp           <- cbind(Xnt_temp, diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])) )
            Xnt_temp <- cbind(Xnt_temp, diagdiff(alt, norm_alt) %*% t(Xnt[, ord[[i_var]], drop = FALSE]))
            # Xnt_temp           <- cbind(Xnt_temp,diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
            colnames(Xnt_temp) <- c(old_col, i_var)
          }
        }
      }

      ### add all the non-alternative specific variables
      for (i_var in allvars[[2]]) {
        if (i_var != 0) {
          old_col <- colnames(Xnt_temp)
          ### remove one regressor (wrt norm_alt), put covariates with random effects at the front
          if (i_var %in% re) {
            Xnt_temp <- cbind(diagn(alt, norm_alt) * Xnt[, i_var], Xnt_temp)
            colnames(Xnt_temp) <- c(paste0(i_var, "_", alt_names), old_col)
          } else {
            Xnt_temp <- cbind(Xnt_temp, diagn(alt, norm_alt) * Xnt[, i_var])
            colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", alt_names))
          }
        }
      }

      ### add all the alternative specific variables with alternative specific coefficient


      for (i_var in allvars[[3]]) {
        if (i_var != 0) {
          old_col <- colnames(Xnt_temp)
          ### put covariates with random effects at the front
          if (i_var %in% re) {
            Xnt_temp <- cbind(diag(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
            # Xnt_temp           <- cbind(diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])), Xnt_temp)
            colnames(Xnt_temp) <- c(paste0(i_var, "_", choice_levels), old_col)
          } else {
            Xnt_temp <- cbind(Xnt_temp, diag(as.vector(Xnt[, ord[[i_var]]])))
            # Xnt_temp           <- cbind(Xnt_temp, diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
            colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", choice_levels))
          }
        }
      }

      ### add alternative specific constants
      if (ASC) {
        old_col <- colnames(Xnt_temp)
        if ("ASC" %in% re) {
          Xnt_temp <- cbind(diagn(alt, norm_alt), Xnt_temp)
          colnames(Xnt_temp) <- c(paste0("ASC", "_", alt_names), old_col)
        } else {
          Xnt_temp <- cbind(Xnt_temp, diagn(alt, norm_alt))
          colnames(Xnt_temp) <- c(old_col, paste0("ASC", "_", alt_names))
        }
      }

      Xnt_temp <- Xnt_temp[, !colnames(Xnt_temp) == "remove_me"]
      Xn[[t]] <- as.matrix(Xnt_temp)
    }

    data[[n]] <- list("X" = Xn, "y" = yn)
  }

  data <- data_cl$new(
    data = data,
    vars = colnames(Xnt_temp),
    ordered = data_raw$ordered
  )

  return(data)
}