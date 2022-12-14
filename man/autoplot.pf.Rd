% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{autoplot.pf}
\alias{autoplot.pf}
\title{Make an automatic \code{ggplot2} plot for a \code{pf} object}
\usage{
\method{autoplot}{pf}(
  object,
  columns = NULL,
  layout = "circular",
  suppress_tiplabels = FALSE,
  ...
)
}
\arguments{
\item{object}{A \code{pf} object to plot}

\item{columns}{Columns to plot along with the phylogeny.
Can use bare column names or any other \code{tidyselect} syntax}

\item{layout}{\code{ggtree::ggtree()} layout to use.}

\item{suppress_tiplabels}{If \code{TRUE}, don't draw tip labels.}

\item{...}{Other arguments passed to or from other methods.}
}
\value{
A \code{ggplot} object.
}
\description{
Make an automatic \code{ggplot2} plot for a \code{pf} object
}
\examples{
autoplot(pf(rpfc(100)) \%>\% dplyr::mutate(trait = rnorm(dplyr::n())), trait,
layout = "rectangular")
}
