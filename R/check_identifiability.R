#' Test for identifiability of a model
#' @description Function that tests whether a given model specification might be identifiable. 
#' @param mod \code{mod}-object
#' @param tol 
#' A real; describes the tolerance in the check. If the determinant changes less than 'tol', the system might be not identifiable. 
#' By default  \code{tol = 0.001}
#' @return boolean ok: FALSE if a reason for non-identifiability is found. 
#' @export
check_identifiability <- function(mod, tol = 0.001){
  ok = FALSE
  npar = mod$lthb + mod$lthL + mod$lthO
  theta = rnorm(npar)
  J = mod$alt
  lthbb = dim(mod$Hb)[1]
  lRE= mod$lRE 
  
  dims = lthbb + lthbb^2 + J^2
  systems = matrix(0,dims,npar)
  
  del1 = diag(J)
  del1[,1]=-1
  sys0 = build_par_from_mod(theta, mod)
  Sig = sys0$Sigma
  b = sys0$b
  Om = sys0$Omega
  
  Sig = del1 %*% Sig %*% t(del1)
  nb = norm(b)
  b = b/norm(b)
  Om = Om/nb^2
  Sig = Sig/nb^2
  vec0 = rbind(b, matrix(Om,ncol=1),matrix(Sig,ncol=1))
  
  for (j in 1:npar){
    thetad = theta
    thetad[j] = thetad[j]+ tol # change each of the coordinates to see, if system changes.
    sysa = build_par_from_mod(thetad, mod)
    Sig = sysa$Sigma
    b = sysa$b
    Om = sysa$Omega
    
    Sig = del1 %*% Sig %*% t(del1)
    nb = norm(b)
    b = b/norm(b)
    Om = Om/nb^2
    Sig = Sig/nb^2
    veca = rbind(b, matrix(Om,ncol=1),matrix(Sig,ncol=1))
    systems[,j] = (veca-vec0)/tol
  }
  
  if (det(t(systems) %*% systems)>0.0001^length(theta)){ok = TRUE}
  
  if (ok == FALSE){
    # check, if scale is not fixed. 
    thetad = theta * (1+tol) # change scale to see, if system changes.
    sysa = build_par_from_mod(thetad, mod)
    Sig = sysa$Sigma
    b = sysa$b
    Om = sysa$Omega
    
    Sig = del1 %*% Sig %*% t(del1)
    nb = norm(b)
    b = b/norm(b)
    Om = Om/nb^2
    Sig = Sig/nb^2
    veca = rbind(b, matrix(Om,ncol=1),matrix(Sig,ncol=1))
    if (norm(vec0-veca)/tol<0.001){
      warning("check_identifiability: Scale does not change system. Consider fixing one entry of beta or Sigma.")
    }
  }                             
  return( ok)
}
