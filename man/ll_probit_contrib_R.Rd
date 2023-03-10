% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ll_probit_contrib_R}
\alias{ll_probit_contrib_R}
\title{Calculation of contribution of one person to the approximated probit log likelihood calculation.}
\usage{
ll_probit_contrib_R(X_n, y_n, theta, mod, approx_method)
}
\arguments{
\item{X_n}{X_n list of regressor matrices}

\item{y_n}{vector of choices}

\item{theta}{parameter vector}

\item{mod}{A \code{\link{mod_cl}} object.}

\item{approx_method}{string specifying the approximation method}
}
\value{
A vector, containing the negative log-likelihood with attributes containing the gradient and (if specified in the controls) the Hessian.
}
\description{
Computes the contribution of one person to the log likelihood.
}
\keyword{internal}
