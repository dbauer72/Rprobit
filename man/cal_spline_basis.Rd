% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cal_spline_basis.R
\name{cal_spline_basis}
\alias{cal_spline_basis}
\title{calculate spline function}
\usage{
cal_spline_basis(grid_points, knots)
}
\arguments{
\item{grid_points}{A real matrix. Each column corresponds to a grid point.}

\item{knots}{A real vector of knots for the spline.}
}
\value{
A real matrix. Each column corresponds to a separate grid point.
}
\description{
This function evaluates the natural cubic splines using knots at the points provided.
}
\keyword{internal}
