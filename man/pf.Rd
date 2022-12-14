% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyf-vctrs.R
\name{pf}
\alias{pf}
\title{\code{pf} object constructor}
\usage{
pf(x = pfc(), pf_column = "phlo", ...)
}
\arguments{
\item{x}{a \code{pfc} object}

\item{pf_column}{Name of the column to hold \code{x}. Default: \code{'phlo'}}

\item{...}{Reserved for future extensions. Not used.}
}
\value{
a \code{pf} object
}
\description{
\code{pf} object constructor
}
\examples{
pf(rpfc(100))
}
