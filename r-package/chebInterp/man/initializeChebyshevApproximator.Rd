% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initializeChebyshevApproximator.R
\name{initializeChebyshevApproximator}
\alias{initializeChebyshevApproximator}
\title{Initializes the Chebyshev Approximation}
\usage{
initializeChebyshevApproximator(D, N, M = N + 1, bounds = NULL,
  upper_b = NULL, lower_b = NULL)
}
\arguments{
\item{D}{Dimensions of the Problem (integer)}

\item{N}{Highest Degree of the Polynomial (integer)}

\item{M}{Number of Interpolation Nodes in each dimension (integer)}

\item{bounds}{Bounds of the rectangle on which the function is approximated (list)}

\item{upper_b}{A vector of upper bounds (numeric)}

\item{lower_b}{A vector of lower bounds (numeric)}
}
\value{
A list of the initialized approximation
}
\description{
Initializes the Chebyshev Approximation
}
