% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_asb.R
\name{plot.asb}
\alias{plot.asb}
\title{Plot asb object}
\usage{
\method{plot}{asb}(
  x,
  parent = "Order",
  child = "Family",
  type = "pa",
  remove_empty_child = FALSE,
  trans = NULL,
  ...
)
}
\arguments{
\item{x}{an object of class asb generatated with \code{as_biomonitor}.}

\item{parent}{name of the coarser taxonomic level.}

\item{child}{name of the finer taxonomic level.}

\item{type}{plot present-absence data with \code{pa}, abundance data with \code{abundance} and
detection probabilities with \code{frequency}.}

\item{remove_empty_child}{remove empty level of the child taxonomic level.}

\item{trans}{transformation for abundance data.}

\item{...}{further arguments to be passed to \code{plot}.}
}
\description{
Plot 2 taxonomic levels using a sunburst visualization.
}
