% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_nonpara_grid.R
\name{fit_nonpara_grid}
\alias{fit_nonpara_grid}
\title{Fit a Probit Model using non-parametric mixing distribution based on a fixed grid.
The parameter for the first regressor is mixed non-parametrically, the remaining parameters are fixed and contained in
the grid elements.}
\usage{
fit_nonpara_grid(data_tr, mod, control, cml_pair_type = 1)
}
\arguments{
\item{data_tr}{data set already preprocessed}

\item{mod}{\code{\link{mod_cl}} object}

\item{control}{\code{\link{control_cl}} object}

\item{cml_pair_type}{An integer. By default, \code{cml_pair_type = 1}. -1 triggers independence likelihood, 0 for FP, 1 for AP.}
}
\value{
A \code{\link{Rprobit_cl}} object.
}
\description{
This function fits a mixed probit model using a non-parametric mixing with fixed grid and estimated mixing
probabilities.
DOES NOT USE PARALLELIZATION!!!
}
