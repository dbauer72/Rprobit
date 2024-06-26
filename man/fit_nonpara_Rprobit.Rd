% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_nonpara_Rprobit.R
\name{fit_nonpara_Rprobit}
\alias{fit_nonpara_Rprobit}
\title{Fit a Probit Model using non-parametric mixing distribution based on a fixed grid.}
\usage{
fit_nonpara_Rprobit(
  Rprobit_obj,
  init_method = "random",
  control_nlm = NULL,
  cml_pair_type = 1,
  iter.adapt = 1
)
}
\arguments{
\item{Rprobit_obj}{An \code{\link{Rprobit_cl}} object.}

\item{init_method}{A \code{character}, setting the method for initialization of the likelihood
optimization.
Can be one of
\itemize{
\item \code{"random"} (default)
\item or \code{"subset"}.
}}

\item{control_nlm}{Optionally a named \code{list} of control parameters passed on to
\link[stats]{nlm}.
By default, \code{control_nlm = NULL}.}

\item{cml_pair_type}{An integer.}

\item{iter.adapt}{an integer specifying the number of adaptations. iter.adapt (default) implies none.
cml_pair_type = 0 (full pairwise), = 1 (adjacent pairwise looped), = -1 (single choices)
By default, \code{cml_pair_type = 1}.}
}
\value{
A \code{\link{Rprobit_cl}} object.
}
\description{
This function fits a mixed probit model using a non-parametric mixing with fixed grid and estimated mixing
probabilities.
DOES NOT USE PARALLELIZATION!!!
}
