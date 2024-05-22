#' Redraw data from Rprobit_cl object
#'
#' @description
#' This function simulates new choices from data stored in an Rprobit_cl object, which contains full information on regressors and model and parameters.
#'
#' @param Rprobit_obj 
#' An \code{\link{Rprobit_cl}} object.
#' @param seed
#' Seed for the simulation of choice data
#' @param thetas 
#' real matrix, each row corresponds to one value of the parameters. When number of persons matches, each person gets one parameter. Else parameter vectors are drawn randomly. 
#'
#' @return list an updated \code{\link{Rprobit_cl}} object. 
#' 

redraw_choices_from_data <- function( Rprobit_obj, seed = 1, thetas = NULL) {
  set.seed(seed)
  
  # decide, if thetas drawn randomly or given
  N <- Rprobit_obj$data$N
  Nth <- dim(thetas)[1]
  rand_theta <- TRUE 

  if (N == Nth){
    rand_theta <- FALSE 
  }
  # clone object not to mix up with input
  Rprobit_o <- Rprobit_obj$clone(deep = TRUE)
  
  # check, if the object is right
  if(!inherits(Rprobit_o, "Rprobit_cl")) {
    stop("'Rprobit_o' is not a Rprobit_cl object")
  }
  if(is.null(Rprobit_o$mod)) {
    stop("'Rprobit_o' does not contain a 'mod'  object")
  }
  
  # check if both data and data_raw are contained 
  if (is.null(Rprobit_o$data) | is.null(Rprobit_o$data_raw)){
    stop("'Rprobit_o' must contain both raw data and adapted data!")
  }
  # switch depending on ordered model or not 
  if (Rprobit_o$mod$ordered == TRUE){
      # ordered data 
      # data_raw contains data frame df; data contains list with (y,X) pairs. 
      # loop over list in data 
      cur <- 0
      for (n in 1:length(Rprobit_o$data$data)){
        
        # choose theta parameter. 
        choice_theta <- n
        if (rand_theta == TRUE){
          choice_theta <- sample(Nth,1)
        }
        # obtain model 
        par <- build_par_from_mod(thetas[choice_theta,], Rprobit_o$mod)
        
        # extract system 
        b <- par$b
        L_Om <- par$Omega_chol
        L_Sig <- par$Sigma_chol
        tauk<- par$tauk
      
        alt <- dim(L_Sig)[1]

        X <- Rprobit_o$data$data[[n]]$X
        y <- matrix(0,dim(X)[1],1)
        Xb <- X %*% b
        utils = Xb + L_Sig %*% stats::rnorm(alt)
        # find categories 
        for (j in 1:dim(X)[1]){
          y[j] <- sum(utils[j] > tauk) + 1
        }
        
        # write data back to object
        Rprobit_o$data_raw$df[cur + c(1:dim(X)[1]),Rprobit_o$data_raw$choice] <- y
        Rprobit_o$data$data[[n]]$y <- y
        cur <- cur + dim(X)[1]
      }
      
  } else {
    # categorical data: data_raw contains data frame, data contains a list of (y,X) data, which itself is a list, one entry for each choice. 
    
    # find levels 
    cats <- unique(Rprobit_o$data_raw$df[,Rprobit_o$data_raw$choice])
    if (is.numeric(cats)){
      lev <- levels(as.factor(cats))
    } else {
      lev <- levels(cats)
    }

    cur <- 0
    for (n in 1:length(Rprobit_o$data$data)){
      
      # choose theta parameter. 
      choice_theta <- n
      if (rand_theta == TRUE){
        choice_theta <- sample(Nth,1)
      }
      # obtain model 
      par <- build_par_from_mod(thetas[choice_theta,], Rprobit_o$mod)
      
      # extract system 
      b <- par$b
      L_Om <- par$Omega_chol
      L_Sig <- par$Sigma_chol
      
      alt <- dim(L_Sig)[1]


      # for each decider extract list of data
      # draw random effects 
      data_list <- Rprobit_o$data$data[[n]]
      y <- data_list$y
      for (j in 1:length(data_list$X)){
        X <- data_list$X[[j]]
        Xb = X %*% b
        utils = Xb + L_Sig %*% stats::rnorm(alt)
        # find best alternative 
        y[[j]] <- which.max(utils)
      }
      
      # write data back to object
      Rprobit_o$data_raw$df[cur + c(1:length(y)),Rprobit_o$data_raw$choice] <- y
      Rprobit_o$data$data[[n]]$y <- y
      cur <- cur + length(y)
    }
  }
  return (Rprobit_o) # return object, although not necessary due to R6 class. 
}

  
