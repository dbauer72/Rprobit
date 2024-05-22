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

data_raw_to_data <- function(data_raw, allvars, choice, re, norm_alt, alt=0, ASC) {
  
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
  choice_levels <- levels(as.factor(data_raw$alt_names))
  
  if (is.null(alt) == TRUE){
    alt = 0
  }
  if (alt == 0){
    alt = length(choice_levels)
  }
  #choice_levels <- data_raw$alt_names
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
  #if(is.numeric(data_raw[,choice])){
  #  alt_names <- levels(as.factor(1:alt))
  # } else {
  #  alt_names <- choice_levels #levels(data_raw[,choice])
  #}
  alt_names <- choice_levels
  
  ordered = data_raw$ordered
  ####################################
  ### split between unordered and ordered data sets, as in the latter case the data recoding works differently. 
  ####################################
  if (ordered == FALSE){
    #### first unordered case ####
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
    
    ### find out the variable names and positions. 
    # this is done for the first decider only 
    n = 1
    
    data_raw_n <- data_raw$df[id_macml == ids[n], ]
    Tp[n] <- dim(data_raw_n)[1]
    yn <- data_raw_n[, choice]
    Xn <- list()
    t = 1
    var_pos = 0
    
    #for (t in 1) {
    Xnt <- data_raw_n[t, ]
    Xnt_temp <- data.frame(remove_me = 1:alt)
    
    DD <- diagdiff(alt, norm_alt)
    DN <- diagn(alt, norm_alt)
    ### add all the alternative specific variables
    for (i_var in allvars[[1]]) {
      if (i_var != 0) {
        old_col <- colnames(Xnt_temp)
        lo = length(ord[[i_var]])
        ### locate the relevant columns and compute difference w.r.t norm_alt, put covariates with random effects at the end
        if (i_var %in% re) {
          # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
          Xnt_temp <-  cbind(DD %*% t(Xnt[, ord[[i_var]], drop = FALSE]), Xnt_temp)
          # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var,"_"), colnames(Xnt))])), Xnt_temp)
          
          
          colnames(Xnt_temp) <- c(i_var, old_col)
          var_pos = c(var_pos+1,1)
        } else {
          # Xnt_temp           <- cbind(Xnt_temp, diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])) )
          Xnt_temp <- cbind(Xnt_temp, DD %*% t(Xnt[, ord[[i_var]], drop = FALSE]))
          # Xnt_temp           <- cbind(Xnt_temp,diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
          colnames(Xnt_temp) <- c(old_col, i_var)
          max_var_pos <- if(length(var_pos)==0) 0 else max(var_pos)
          var_pos = c(var_pos,max_var_pos+1)
        }
      }
    }
    
    ### add all the non-alternative specific variables
    for (i_var in allvars[[2]]) {
      if (i_var != 0) {
        old_col <- colnames(Xnt_temp)
        lo = length(alt_names)
        ### remove one regressor (wrt norm_alt), put covariates with random effects at the front
        if (i_var %in% re) {
          Xnt_temp <- cbind(DN * Xnt[, i_var], Xnt_temp)
          colnames(Xnt_temp) <- c(paste0(i_var, "_", alt_names), old_col)
          var_pos = c(var_pos+lo,1:lo)
        } else {
          Xnt_temp <- cbind(Xnt_temp, DN * Xnt[, i_var])
          colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", alt_names))
          var_pos = c(var_pos,max(var_pos)+c(1:lo))
        }
      }
    }
    
    ### add all the alternative specific variables with alternative specific coefficient
    
    
    for (i_var in allvars[[3]]) {
      if (i_var != 0) {
        old_col <- colnames(Xnt_temp)
        lo = length(ord[[i_var]])
        ### put covariates with random effects at the front
        if (i_var %in% re) {
          Xnt_temp <- cbind(diag(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
          # Xnt_temp           <- cbind(diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])), Xnt_temp)
          colnames(Xnt_temp) <- c(paste0(i_var, "_", choice_levels), old_col)
          var_pos = c(var_pos+lo,1:lo)
        } else {
          Xnt_temp <- cbind(Xnt_temp, diag(as.vector(Xnt[, ord[[i_var]]])))
          # Xnt_temp           <- cbind(Xnt_temp, diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
          colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", choice_levels))
          var_pos = c(var_pos,max(var_pos)+c(1:lo))
        }
      }
    }
    
    ### add alternative specific constants
    if (ASC) {
      old_col <- colnames(Xnt_temp)
      lo = length(alt_names)
      if ("ASC" %in% re) {
        Xnt_temp <- cbind(DN, Xnt_temp)
        colnames(Xnt_temp) <- c(paste0("ASC", "_", alt_names), old_col)
        var_pos = c(var_pos+lo,1:lo)
      } else {
        Xnt_temp <- cbind(Xnt_temp,DN)
        colnames(Xnt_temp) <- c(old_col, paste0("ASC", "_", alt_names))
        var_pos = c(var_pos,max(var_pos)+c(1:lo))
      }
    }
    
    Xnt_temp <- Xnt_temp[, !colnames(Xnt_temp) == "remove_me",drop=FALSE]
    col_names_tot <- colnames(Xnt_temp)
    
    var_pos = var_pos[-1]
    
    ### loop over decision makers
    for (n in 1:N) {
      data_raw_n <- data_raw$df[id_macml == ids[n], ]
      Tp[n] <- dim(data_raw_n)[1]
      yn <- data_raw_n[, choice]
      Xn <- list(Tp[n])
      
      ### loop over choice occasions
      for (t in 1:Tp[n]) {
        Xnt <- data_raw_n[t, ]
        
        Xnt_temp <- matrix(NA,alt,length(col_names_tot))
        #colnames(Xnt_temp) <- col_names_tot
        cur = 1
        ### add all the alternative specific variables
        for (i_var in allvars[[1]]) {
          if (i_var != 0) {
            #old_col <- colnames(Xnt_temp)
            ### locate the relevant columns and compute difference w.r.t norm_alt, put covariates with random effects at the end
            #if (i_var %in% re) {
            # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
            # Xnt_temp <- cbind(diagdiff(alt, norm_alt) %*% t(Xnt[, ord[[i_var]], drop = FALSE]), Xnt_temp)
            Xnt_temp[,var_pos[cur]] <- DD %*% t(Xnt[, ord[[i_var]], drop = FALSE])
            cur = cur+1
            # Xnt_temp           <- cbind(diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var,"_"), colnames(Xnt))])), Xnt_temp)
            #  colnames(Xnt_temp) <- c(i_var, old_col)
            #} else {
            # Xnt_temp           <- cbind(Xnt_temp, diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, ord[[i_var]]])) )
            #  Xnt_temp <- cbind(Xnt_temp, diagdiff(alt, norm_alt) %*% t(Xnt[, ord[[i_var]], drop = FALSE]))
            # Xnt_temp           <- cbind(Xnt_temp,diagdiff(alt, norm_alt)%*%t(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
            #  colnames(Xnt_temp) <- c(old_col, i_var)
            #}
          }
        }
        
        ### add all the non-alternative specific variables
        for (i_var in allvars[[2]]) {
          if (i_var != 0) {
            #old_col <- colnames(Xnt_temp)
            ### remove one regressor (wrt norm_alt), put covariates with random effects at the front
            #if (i_var %in% re) {
            lo = length(alt_names)
            Xnt_temp[,var_pos[cur+ c(0:(lo-1))]] <- DN * Xnt[, i_var] #cbind(diagn(alt, norm_alt) * Xnt[, i_var], Xnt_temp)
            cur = cur+lo
            #  colnames(Xnt_temp) <- c(paste0(i_var, "_", alt_names), old_col)
            #} else {
            #  Xnt_temp <- cbind(Xnt_temp, diagn(alt, norm_alt) * Xnt[, i_var])
            #  colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", alt_names))
            #}
          }
        }
        
        ### add all the alternative specific variables with alternative specific coefficient
        
        
        for (i_var in allvars[[3]]) {
          if (i_var != 0) {
            #old_col <- colnames(Xnt_temp)
            lo = length(choice_levels)
            ### put covariates with random effects at the front
            #if (i_var %in% re) {
            Xnt_temp[,var_pos[cur+ c(0:(lo-1))]] <- diag(as.vector(Xnt[, ord[[i_var]]]))    #    cbind(diag(as.vector(Xnt[, ord[[i_var]]])), Xnt_temp)
            cur = cur+lo
            # Xnt_temp           <- cbind(diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])), Xnt_temp)
            #  colnames(Xnt_temp) <- c(paste0(i_var, "_", choice_levels), old_col)
            #} else {
            #  Xnt_temp <- cbind(Xnt_temp, diag(as.vector(Xnt[, ord[[i_var]]])))
            #  # Xnt_temp           <- cbind(Xnt_temp, diag(as.vector(Xnt[, grep(paste0(i_var, "_"), colnames(Xnt))])))
            #  colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", choice_levels))
            #}
          }
        }
        
        ### add alternative specific constants
        if (ASC) {
          #old_col <- colnames(Xnt_temp)
          #if ("ASC" %in% re) {
          lo = alt-1
          Xnt_temp[,var_pos[cur+ c(0:(lo-1))]] <- DN  # cbind(diagn(alt, norm_alt), Xnt_temp)
          cur = cur+lo
          #colnames(Xnt_temp) <- c(paste0("ASC", "_", alt_names), old_col)
          #} else {
          #  Xnt_temp <- cbind(Xnt_temp, diagn(alt, norm_alt))
          #  colnames(Xnt_temp) <- c(old_col, paste0("ASC", "_", alt_names))
          #}
        }
        
        #Xnt_temp <- Xnt_temp[, !colnames(Xnt_temp) == "remove_me"]
        Xnt_df = data.frame(Xnt_temp)
        colnames(Xnt_df) <- col_names_tot
        Xn[[t]] <- Xnt_temp
      }
      
      data[[n]] <- list("X" = Xn, "y" = yn)
    }
    
    
    
  } else {
    #################################
    #### now the ordered case 
    #################################
    
    if (!is.numeric(norm_alt)) {
      norm_alt <- which(alt_names == norm_alt)
    }
    
    if (norm_alt != 0) alt_names <- alt_names[-norm_alt]
    data <- list()
    
    data_raw_n <- data_raw$df[data_raw$df[, id] == ids[1], ]
    Xnt <- data_raw_n[1, ]
    
    ord <- list()
    ord[["ASC"]] <- c() # no ASCs in ordered probit. 
    
    
    ### find out the variable names and positions. 
    # this is done for the first decider only 
    n = 1
    
    data_raw_n <- data_raw$df[id_macml == ids[n], ]
    Tp[n] <- dim(data_raw_n)[1]
    yn <- data_raw_n[, choice]
    Xn <- list()
    var_pos = 0
    
    #for (t in 1) {
    Xnt <- data_raw_n[1,]
    Xnt_temp <- data.frame(remove_me = 1)
    
    ### add all the non-alternative specific variables
    for (i_var in allvars[[2]]) {
      if (i_var != 0) {
        old_col <- colnames(Xnt_temp)
        lo = length(alt_names)
        ### remove one regressor (wrt norm_alt), put covariates with random effects at the front
        if (i_var %in% re) {
          Xnt_temp <- cbind(Xnt[1, i_var], Xnt_temp)
          colnames(Xnt_temp) <- c(i_var, old_col)
          var_pos = c(var_pos+1,1)
        } else {
          Xnt_temp <- cbind(Xnt_temp, Xnt[1, i_var])
          colnames(Xnt_temp) <- c(old_col, i_var)
          var_pos = c(var_pos,max(var_pos)+1)
        }
      }
    }
    
    Xnt_temp <- Xnt_temp[, !colnames(Xnt_temp) == "remove_me"]
    col_names_tot <- colnames(Xnt_temp)
    
    var_pos = var_pos[-1]
    
    
    ### loop over decision makers
    for (n in 1:N) {
      data_raw_n <- data_raw$df[id_macml == ids[n], ]
      Tp <- dim(data_raw_n)[1]
      yn <- data_raw_n[, choice]
      Xnt <- data_raw_n
      
      Xnt_temp <- matrix(NA,Tp,dim(Xnt_temp)[2])
      colnames(Xnt_temp) <- col_names_tot
      
      ### add all the non-alternative specific variables
      cur = 1
      for (i_var in allvars[[2]]) {
        if (i_var != 0) {
          #old_col <- colnames(Xnt_temp)
          ### remove one regressor (wrt norm_alt), put covariates with random effects at the front
          #if (i_var %in% re) {
          Xnt_temp[,var_pos[cur]] <- Xnt[, i_var] #cbind(diagn(alt, norm_alt) * Xnt[, i_var], Xnt_temp)
          cur = cur+1
          #  colnames(Xnt_temp) <- c(paste0(i_var, "_", alt_names), old_col)
          #} else {
          #  Xnt_temp <- cbind(Xnt_temp, diagn(alt, norm_alt) * Xnt[, i_var])
          #  colnames(Xnt_temp) <- c(old_col, paste0(i_var, "_", alt_names))
          #}
        }
      }
      
      Xn <- Xnt_temp
      Xnt_df = data.frame(Xnt_temp)
      colnames(Xnt_df) <- col_names_tot
      
      if (inherits(data_raw,"data_raw_StSp_cl")){
        data[[n]] <- list("X" = Xn, "y" = yn, "time" = data_raw_n[,data_raw$time_col], "quest" = data_raw_n[,data_raw$quest_col])
      } else {
        data[[n]] <- list("X" = Xn, "y" = yn)
      }     
    } # end cycle over deciders 
  } # end switch for ordered and unordered. 
  if (inherits(data_raw,"data_raw_StSp_cl")){
    data_obj <- data_StSp_cl$new(
      data = data,
      vars = colnames(Xnt_df),
      ordered = data_raw$ordered
    )
    data_obj$time <- unique(data_raw$df[,data_raw$time_col])
    data_obj$quest = max(data_raw$df[,data_raw$quest_col])
    
  } else {
    data_obj <- data_cl$new(
      data = data,
      vars = colnames(Xnt_df),
      ordered = data_raw$ordered
    )
  }
  
  return(data_obj)
}


