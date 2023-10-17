#' Plot the coefficient distribution function 
#' (makes sense, if random coefficients or latent classes exist) 
#'
#' @description
#' This function provides a plot for the 
#'
#' @inheritParams AIC.Rprobit_cl
#' @param dims 
#' vector of integers indicating the coordinates to be plotted. Only the first two integers are considered
#' @param cdf 
#' Boolean, indicating whether the cdf (TRUE) or the pdf (FALSE; default) should be plotted.  
#' @param margin
#' real; factor (0.1 = 10%; minimum default margin) to enlarge plotting area. 
#' @return
#' A plot.
#'
#' @export

plot_Rprobit <- function(object, dims = 1, cdf = FALSE, margin = 0.1) {
  # check, if object is latent class or contains random coefficients. 
  LC_model = FALSE
  ### build model parameters
  if (!is.null(object$mod$num_class)){
    LC_model = TRUE
    num_class <- object$mod$num_class 
  } else {
    num_class <- 1
  }
  
  lRE <- object$mod$lRE 
  
  if ((LC_model == FALSE) && (lRE == 0)){
    paste("The model does not contain random effects and it is not a latent class model. The plot will be really boring!")
  }
  
  par_all <- list(num_class)
  mod <- object$mod
  for (j in 1:num_class){
    param_one <- mod$lthb + mod$lthO + mod$lthL 
    ind_j <- (j-1)*param_one+c(1:param_one)
    par_all[[j]] <- build_par_from_mod(object$theta[ind_j], mod)
  }
  if (num_class>1){ 
    pi_eta <- c(0,object$theta[num_class*param_one + c(1:(num_class-1))])
    pi <-  exp(pi_eta)/sum(exp(pi_eta))
  } else {
    pi <- 1
  }
  
  # one dimensional plot
  if (length(dims) == 1){
    # compile data 
    beta <- matrix(0,num_class,1)
    omega <- beta
    for (j in 1:num_class){
      beta[j] <- par_all[[j]]$b[dims]
      omega[j] <- par_all[[j]]$Omega[dims,dims]
    }
    sd_om = sqrt(omega)
    range_beta <- c(min(beta),max(beta))
    drb <- diff(range_beta)
    if (drb<0.1){ drb <- 0.1}
    if (margin<0.1){margin <- 0.1}

    range_beta <- c(range_beta[1]- margin*drb , range_beta[2]+ margin*drb) 
    xgrid <- seq(from=range_beta[1],to=range_beta[2],by= diff(range_beta)/200)
    if (max(omega)>0){
      if (cdf == TRUE){
        curv <- stats::pnorm(xgrid,mean = beta[1],sd = sd_om[1])*pi[1]
        curv_tot <- curv
        graphics::plot(xgrid,curv,col=1,type='l',ylim = c(0,1))
        if (num_class>1){
          for (j in 2:num_class){
            curv <- stats::pnorm(xgrid,mean = beta[j],sd = sd_om[j])*pi[j]
            curv_tot <- curv_tot+  curv
            graphics::lines(xgrid,curv,col=j)
          }
        }
        graphics::lines(xgrid,curv_tot,col="black",cex = 2)
      } else {
        # pdf wanted. 
        curv <- stats::dnorm(xgrid,mean = beta[1],sd = sd_om[1])*pi[1]
        curv_tot <- curv
        ymax <- max(stats::dnorm(0,mean=0,sd = min(sd_om)))
        graphics::plot(xgrid,curv,col=1,type='l',ylim = c(0,ymax))
        if (num_class>1){
          for (j in 2:num_class){
            curv <- stats::dnorm(xgrid,mean = beta[j],sd = sd_om[j])*pi[j]
            curv_tot <- curv_tot+  curv
            graphics::lines(xgrid,curv,col=j)
          }
        }
        graphics::lines(xgrid,curv_tot,col="black",cex = 2)
      }
    } else {
      # no random effects -> different plots 
      if (cdf == TRUE){
        ord_beta <- base::order(beta)
        sfun0  <- stats::stepfun(beta[ord_beta], c(0,cumsum(pi[ord_beta])), f = 0)
        graphics::plot(sfun0, main = "CDF of coefficient")
      } else {
        # PDF wanted
        graphics::plot(beta,pi)
      }
    }
    
    # find ranges 
    
    
  } else {
    if (length(dims)>2){
      paste("Only the first two components are plotted.")
    }
    dims = dims[c(1,2)]
    # get data 
    beta <- matrix(0,num_class,2)
    omega <- list(num_class)
    sd_omega <- beta
    for (j in 1:num_class){
      beta[j,] <- par_all[[j]]$b[dims]
      omega[[j]] <- par_all[[j]]$Omega[dims,dims]
      sd_omega[j,] <- sqrt(diag(omega[[j]]))
    }
    range_beta_x <- c(min(beta[,1]),max(beta[,1]))
    drbx <- diff(range_beta_x)
    if (drbx<0.1){drbx <- 0.1}
    if (margin<0.1){margin <- 0.1}
    range_beta_x <- c(range_beta_x[1]- drbx*margin , range_beta_x[2]+ drbx*margin) 
    
    range_beta_y <- c(min(beta[,2]),max(beta[,2]))
    drby <- diff(range_beta_y)
    if (drby<0.1){drby <- 0.1}
    range_beta_y <- c(range_beta_y[1]- drby*margin , range_beta_y[2]+ drby*margin) 
    
    xgrid <- seq(from=range_beta_x[1],to=range_beta_x[2],by= diff(range_beta_x)/50)
    ygrid <- seq(from=range_beta_y[1],to=range_beta_y[2],by= diff(range_beta_y)/50)
    if (max(sd_omega)>0){
      # there is at least one random component
      sigma <- matrix(c(1, 0, 0, 1), nrow=2)
      if (cdf == TRUE){
        f     <- function(x, y) {
          fn <- matrix(0,length(x),1)
          for (j in 1:num_class){
            om <- omega[[j]]
            if (om[1,1]<0.0001){ om[1,1] = 0.0001}
            if (om[2,2]<0.0001){ om[2,2] = 0.0001}
            rho = om[1,2]/sqrt(om[1,1]*om[2,2])
            som1 = sqrt(om[1,1])
            som2 = sqrt(om[2,2])
            for (n in 1:length(x)){
              fn[n] <- fn[n] + pi[j]*biv_normal_cdf((x[n]-beta[j,1])/som1,(y[n]-beta[j,2])/som2,rho)
            }
          }
          return (fn)
        }
      } else {
        f     <- function(x, y) {
          fn <- matrix(0,length(x),1)
          for (j in 1:num_class){
            om <- omega[[j]]
            if (om[1,1]<0.0001){ om[1,1] = 0.0001}
            if (om[2,2]<0.0001){ om[2,2] = 0.0001}
            rho = om[1,2]/sqrt(om[1,1]*om[2,2])
            som1 = sqrt(om[1,1])
            som2 = sqrt(om[2,2])
            for (n in 1:length(x)){
              fn[n] <- fn[n] + pi[j]*biv_normal_pdf((x[n]-beta[j,1])/som1,(y[n]-beta[j,2])/som2,rho)
            }
          }
          return (fn)
        }
      }
      z     <- base::outer(xgrid, ygrid, f)
    
      #create contour plot
      graphics::contour(xgrid, ygrid, z, add = FALSE, col = "black")
    } else {
      graphics::plot(beta[,1],beta[,2],xlim = range_beta_x,ylim = range_beta_y)
    }
      
    # add means and percentages to plot 

    pi_norm = pi/max(pi)
    for (j in 1:num_class){
      graphics::points(beta[j,1],beta[j,2],cex=2*pi_norm[j], pch = 19,col=j)
      graphics::text(beta[j,1],beta[j,2], sprintf("[%1.2f]",pi[j]),cex=1.2, pos=3,col=j)
    }
  }

  
} 
