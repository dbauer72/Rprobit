#' Simulate data_raw
#'
#' @description
#' This function simulates choice data.
#'
#' @param ASC
#' Boolean, indicating whether ASCs are estimated or not
#' @param allvars
#' \code{allvars}-object
#' @param re
#' Character vector of variable names which are considered to have random effects
#' @param choice
#' Character, name of dependent variable
#' @param theta
#' \code{theta}-object
#' @param mod
#' A \code{\link{mod_cl}} or \code{\link{mod_latclass_cl}} object.
#' @param seed
#' Seed for the simulation of choice data
#' @param control_simulation
#' list to fine-tune the simulated data. In case auto correlated errors should be simulated, control_simulation$simulate_ar1_rho determines auto correlation coefficient/matrix
#'
#' @return list of \code{data_raw}, \code{data} and \code{theta}
#' 
#' @keywords internal

sim_data_raw <- function(
    ASC, allvars, re, choice, theta, mod, seed, control_simulation = NULL
  ) {

  ### differentiate between latent class and else 
  lat_class = FALSE
  num_class = 1 
  
  if (inherits(mod,"mod_latclass_cl")==TRUE){
    lat_class = TRUE
    num_class = mod$num_class
  }
  ### set seed
  if (!is.null(seed)) set.seed(seed)

  ### draw parameter vector
  if (is.null(theta)) {
    theta <- stats::rnorm(mod$lthb + mod$lthO + mod$lthL)

    if (mod$ordered == TRUE) {
      dtau <- stats::rnorm(mod$alt - 1, 0, 1)
      tauk <- (exp(dtau))
      theta <- c(theta, tauk)
    }
    if (lat_class==TRUE){
      # latent class model requires more parameters. 
      tot_params <- mod$tot_params
      theta <- stats::rnorm(tot_params)
    }
  }

  ### build 'par'
  # list in case latent classes are wanted. 
  par_all <- list(num_class)
  if (lat_class == FALSE){
    par_all[[1]] <- build_par_from_mod(theta, mod)
    pi = 1
  } else {
    for (j in 1:num_class){
      param_one <- mod$lthb + mod$lthO + mod$lthL 
      par_all[[j]] <- build_par_from_mod(theta[(j-1)*param_one+c(1:param_one)], mod)
    }
    pi_eta = c(0,theta[num_class*param_one + c(1:(num_class-1))])
    pi = exp(pi_eta)/sum(exp(pi_eta))
  }
  
  
  
  ### extract parameters
  
  if (!is.null(control_simulation)) {
    Tp <- control_simulation$Tp
    N <- length(Tp) # mod$N
    ### distribution used for random effects
    if (!is.null(control_simulation$dis)){
      dis = control_simulation$dis
    } else {
      dis = 'rnorm('
    }
    if (!is.null(control_simulation$pars)){
      pars = control_simulation$pars
    } else {
      pars = ',0,1)'
    }
    
  } else {
    Tp <- rep(1, 100)
    N <- 100
    ### distribution used for random effects
    dis = 'rnorm('
    pars = ',0,1)'
  }
  alt <- mod$alt

  ### draw class membership
  member_draws = stats::rmultinom(N,1,pi)
  class_member <- as.list((1:num_class) %*% member_draws)
  
  
  ### extract simulation controls
  if (!is.null(control_simulation)) {
    rho <- control_simulation$simulate_ar1_rho
  } else {
    rho <- 0
  }
  ### split between ordered and unordered probit model for data generation
  if (mod$ordered == TRUE) {
    if ((N == 1) || (stats::var(Tp) > 0)) warning("For ordered probit Tp must be identical for all persons. Only first value is respected.")
    Tp <- Tp[1]
    data_raw_df <- data.frame()
    data_raw <- data_raw_cl$new(choice = "choice", id = "id_macml")

    ### extract category bounds
    par <- par_all[[1]]
    tauk <- par$tauk

    for (n in 1:N) {

      par <- par_all[[ class_member[[n]] ]]
      b <- par$b
      Omega_chol <- par$Omega_chol
      Sigma_chol <- par$Sigma_chol
      
      ### linear coefficient vector
      b_n <- b + Omega_chol %*% eval(parse(text = sprintf('stats::%s %d%s',dis,length(b),pars)))    # stats::rnorm(length(b))
      e_n <- Sigma_chol %*% stats::rnorm(Tp)

      ### type 1 covariate
      Xn <- matrix(stats::rnorm(Tp * length(allvars[[2]]), 0, 1), nrow = Tp, ncol = length(allvars[[2]]))
      colnames(Xn) <- allvars[[2]]

      utils <- Xn %*% b_n + e_n
      yn <- matrix(0, Tp, 1)
      for (j in 1:Tp) {
        yn[j] <- sum(utils[j] > tauk) + 1
      }
      data_raw_net <- data.frame(choice = yn)
      data_raw_net <- cbind(data_raw_net, Xn)
      data_raw_net$"id_macml" <- n
      data_raw_df <- rbind(data_raw_df, data_raw_net)
    }
    ### end for n in 1:N
    data_raw$set_df(data_raw_df)
    data <- data_raw_to_data(data_raw = data_raw, allvars = allvars, choice = choice, re = re, norm_alt = 1, alt = alt, ASC = ASC)    
    
  } else {
    ### check 'Tp'
    # if(!length(Tp) %in% c(1, mod$N)) stop("Variable 'Tp' has to be of length 1 or of length equal to number of decision makers.")
    # if(length(Tp)==1) Tp <- rep(Tp, mod$N)

    ### build data_raw
    data_raw_df_li <- list(N) 
    data_raw <- data_raw_cl$new()
    data_raw$id <- "id_macml"
    data_raw$choice <- "choice"
    
    for (n in 1:N) {

      #par <- par_all[[ class_member[[n]] ]]
      #b <- par$b
      #Omega_chol <- par$Omega_chol
      #Sigma_chol <- par$Sigma_chol
      
      ### draw 
      data_raw_nt <- data.frame(id_macml = rep(1,Tp[n]), choice = NA)
      ### linear coefficient vector
      #b_n <- b + Omega_chol %*% stats::rnorm(length(b))

      for (t in 1:Tp[n]) {

        ### characteristics
        cov_nt <- matrix(NA, nrow = 1, ncol = 0)
        ### type 1 covariate
        if (sum(as.character(allvars[[1]]) != "0") > 0) {

          ### check if variable distributions were defined
          if (!is.null(control_simulation$variate_distribution)) {
            allvars_1_names <- allvars[[1]][as.character(allvars[[1]]) != "0"]
            new_data <- c()
            ### go through all variables of type 1
            for (i_vars in 1:length(allvars_1_names)) {
              if (!is.null(control_simulation$variate_distribution[[allvars_1_names[i_vars]]])) {
                new_data_i <- control_simulation$variate_distribution[[allvars_1_names[i_vars]]](alt)
              } else {
                new_data_i <- stats::rnorm(alt)
              }
              new_data <- c(new_data, new_data_i)
            }
          } else {
            new_data <- stats::rnorm(alt * sum(as.character(allvars[[1]]) != "0"))
          }

          cov_nt <- cbind(
            cov_nt,
            matrix(
              data = new_data,
              nrow = 1,
              dimnames = list(1, do.call(paste0, expand.grid(paste0(unlist(allvars[[1]][as.character(allvars[[1]]) != "0"]), "_"), 1:alt)))
            )
          )
        }

        ### type 2 covariates
        # cov_nt <- cbind(cov_nt, matrix(stats::rnorm(length(allvars[[2]])), nrow=1, dimnames = list(1,unlist(allvars[[2]]))))
        if (sum(as.character(allvars[[2]]) != "0") > 0) {

          ### check if variable distributions were defined
          if (!is.null(control_simulation$variate_distribution)) {
            allvars_2_names <- allvars[[2]][as.character(allvars[[2]]) != "0"]
            new_data <- c()
            ### go through all variables of type 1
            for (i_vars in 1:length(allvars_2_names)) {
              if (!is.null(control_simulation$variate_distribution[[allvars_2_names[i_vars]]])) {
                new_data_i <- control_simulation$variate_distribution[[allvars_2_names[i_vars]]](1)
              } else {
                new_data_i <- stats::rnorm(1)
              }
              new_data <- c(new_data, new_data_i)
            }
          } else {
            new_data <- stats::rnorm(sum(as.character(allvars[[2]]) != "0"))
          }

          cov_nt <- cbind(
            cov_nt,
            matrix(
              data = new_data,
              nrow = 1,
              dimnames = list(1, unlist(allvars[[2]][as.character(allvars[[2]]) != "0"]))
            )
          )
        }

        ### type 3 covariate
        # cov_nt <- cbind(cov_nt, matrix(stats::rnorm(alt*length(allvars[[3]])),nrow=1,dimnames=list(1,do.call(paste0,expand.grid(paste0(unlist(allvars[[3]]),"_"), 1:alt)))))
        if (sum(as.character(allvars[[3]]) != "0") > 0) {

          ### check if variable distributions were defined
          if (!is.null(control_simulation$variate_distribution)) {
            allvars_3_names <- allvars[[3]][as.character(allvars[[3]]) != "0"]
            new_data <- c()
            ### go through all variables of type 1
            for (i_vars in 1:length(allvars_3_names)) {
              if (!is.null(control_simulation$variate_distribution[[allvars_3_names[i_vars]]])) {
                new_data_i <- control_simulation$variate_distribution[[allvars_3_names[i_vars]]](alt)
              } else {
                new_data_i <- stats::rnorm(alt)
              }
              new_data <- c(new_data, new_data_i)
            }
          } else {
            new_data <- stats::rnorm(alt * sum(as.character(allvars[[3]]) != "0"))
          }

          cov_nt <- cbind(
            cov_nt,
            matrix(
              data = new_data,
              nrow = 1,
              dimnames = list(1, do.call(paste0, expand.grid(paste0(unlist(allvars[[3]][as.character(allvars[[3]]) != "0"]), "_"), 1:alt)))
            )
          )
        }

        ### create regressor matrix
        data_raw_nt[t,colnames(cov_nt)] <- cov_nt
        #data_raw_nt[t,] <- c(id_macml = 1, choice = NA, cov_nt)
        #colnames(data_raw_nt)[2] <- choice

        
        ### fill data_raw
        data_raw_nt[t,"id_macml"] <- n
        data_raw_nt[t,"choice"] <- 1
        
      }
      
      data_raw_df_li[[n]] <- data_raw_nt
    } ### end for n in 1:N
    
    #data_raw_df <- dplyr::bind_rows(data_raw_df_li)
    data_raw_df <- do.call(rbind, data_raw_df_li)
    data_raw$id <- "id_macml"
    data_raw$choice <- "choice"
    data_raw$set_df(data_raw_df)
         
        #data_raw_df <- rbind(data_raw_df, data_raw_nt)
    # convert raw_data to data for drawing choices. 
    data <- data_raw_to_data(data_raw = data_raw, allvars = allvars, choice = choice, re = re, norm_alt = 1, alt = alt, ASC = ASC)    
        
    ### The following code is there to draw choice from utility. 
    cur = 1
    for (n in 1:N) {
      #data_raw_nt <- data.frame(id_macml = rep(1,Tp[n]), choice = NA)
      par <- par_all[[ class_member[[n]] ]]
      b <- par$b
      Omega_chol <- par$Omega_chol
      Sigma_chol <- par$Sigma_chol
      
      ### linear coefficient vector
      b_n <- b + Omega_chol %*% eval(parse(text = sprintf('stats::%s %d%s',dis,length(b),pars)))  # %*% stats::rnorm(length(b))
      
      for (t in 1:Tp[n]) {
    
        #data_raw$set_df(data_raw_nt[t,])
        #X_nt <- data_raw_to_data(data_raw = data_raw, allvars = allvars, choice = choice, re = re, norm_alt = 1, alt = alt, ASC = ASC)$data[[1]]$X[[1]]
        X_nt <- data$data[[n]]$X[[t]]

        ### draw random error
        error_nt <- Sigma_chol %*% stats::rnorm(alt)
        ### calculate error due to auto correlation
        if (!is.null(rho)) {

          ### in case of simple auto correlation
          if (is.numeric(rho) & is.null(dim(rho))) {
            if (t == 1) {
              error_nt <- error_nt / sqrt(1 - rho^2)
            } else {
              if (!is.null(control_simulation$positions)) {
                positions_n <- control_simulation$positions[[n]]

                t_current <- positions_n[t - 1] + 1
                while (t_current < positions_n[t]) {
                  error_nt <- rho * error_ntm1 + error_nt
                  error_ntm1 <- error_nt
                  error_nt <- Sigma_chol %*% stats::rnorm(alt)
                  t_current <- t_current + 1
                }
              }
              error_nt <- rho * error_ntm1 + error_nt
            }
          } else if (all(dim(rho) == c(max(dim(error_nt)), max(dim(error_nt))))) {
            ### TODO: implement general form for auto correlation
            if (t == 1) {
              kron_rho <- kronecker(rho, rho)
              Sigma_n_tilde <- (Sigma_chol) %*% t(Sigma_chol)
              Sigma_n_tilde_vec <- matrixcalc::vec(Sigma_n_tilde)
              Sigma_n_vec <- solve(diag(nrow(kron_rho)) - kron_rho) %*% Sigma_n_tilde_vec
              Sigma_n <- matrix(Sigma_n_vec, nrow = sqrt(nrow(Sigma_n_vec)))
              Sigma_n_chol <- t(chol(Sigma_n))
              error_nt <- Sigma_n_chol %*% error_nt
            } else {
              if (!is.null(control_simulation$positions)) {
                positions_n <- control_simulation$positions[[n]]

                t_current <- positions_n[t - 1] + 1
                while (t_current < positions_n[t]) {
                  error_nt <- rho %*% error_ntm1 + error_nt
                  error_ntm1 <- error_nt
                  error_nt <- Sigma_chol %*% stats::rnorm(alt)
                  t_current <- t_current + 1
                }
              }
              error_nt <- rho %*% error_ntm1 + error_nt
            }
          }
        }


        ### determine choice
        choice_nt <- which.max(X_nt %*% b_n + error_nt)

        ### save error for next observation for auto regression model
        error_ntm1 <- error_nt

        ### fill data_raw
        data_raw$df[cur,"choice"] <- choice_nt
        data$data[[n]]$y[t] <- choice_nt
        cur = cur+1 
        #data_raw_nt[t,"choice"] <- choice_nt
        #data_raw_df <- rbind(data_raw_df, data_raw_nt)
      }
      
      #data_raw_df_li[[n]] <- data_raw_nt
    } ### end for n in 1:N
    
    #data_raw_df <- bind_rows(data_raw_df_li)
    #data_raw$id <- "id_macml"
    #data_raw$choice <- "choice"
    #data_raw$set_df(data_raw_df)
  }
  data_raw$class_member <- class_member
  data$class_member <- class_member
  
  return(list("data_raw" = data_raw, "data" = data, "theta_0" = theta))
}
