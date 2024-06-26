% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make.path.R
\name{make.path}
\alias{make.path}
\title{Generation of estimated paths of coefficients for Maxent, GM, rGM and Fisher}
\usage{
make.path(
  env = Env2,
  sp = sp_pteridium_aquilinum,
  iter = 500,
  Tau = c(10^{
     -4:2
 })
)
}
\arguments{
\item{sp}{index vector for presence locations.}

\item{iter}{maximum number of maximum iterations.}

\item{Xb}{data matrix of environmental variables and bias variables.}
}
\value{
estimated paths of coefficients
}
\description{
\code{make.path} returns the estimated habitat maps by \code{\link{maxent}}, \code{\link{rgm}}, \code{\link{gm}} and \code{\link{fisher0}}.
}
\examples{
# Example usage
make.path()
}
