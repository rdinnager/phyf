% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyf-vctrs.R
\name{pf_as_pf}
\alias{pf_as_pf}
\title{Convert an object to a \code{pf} object}
\usage{
pf_as_pf(x, ...)
}
\arguments{
\item{x}{An object to convert}

\item{...}{Other arguments pass to or from other methods.}
}
\value{
a \code{pf} object with branch lengths as features
}
\description{
Convert an object to a \code{pf} object
}
\examples{
pf_as_pf(ape::rtree(100))
}
