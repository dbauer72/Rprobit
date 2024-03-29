% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ll_macml_LC}
\alias{ll_macml_LC}
\title{Calculation of Approximated Composite Marginal Likelihood for latent class models. No gradient supplied!}
\usage{
ll_macml_LC(theta, data_obj, mod, control)
}
\arguments{
\item{theta}{parameter vector}

\item{data_obj}{data_cl object}

\item{mod}{A \code{\link{mod_cl}} object.}

\item{control}{Controls for the probit estimation.}
}
\value{
A vector, containing the negative log-likelihood, no gradient and Hessian.
}
\description{
This function computes the composite marginal likelihood for latent class model.
}
