% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manage_traits.R
\name{manage_traits}
\alias{manage_traits}
\title{manage_traits}
\usage{
manage_traits(x, method = "nearest+-", traceB = FALSE)
}
\arguments{
\item{x}{Result of \code{assign_traits()}.}

\item{method}{Can be \code{nearest}, \verb{nearest+}, \verb{nearest-}, \verb{nearest+-} and \verb{neareast-+}.
Please see details for further information.}

\item{traceB}{If \code{TRUE} it will return a vector containing taxa excluded from the selection process
because they did not meet the requirments of the selection.}
}
\description{
A function to select traits based on taxonomic distance.
}
\details{
Method \code{nearest} selects the traits belonging to the nearest taxa irrispective of their position
in the taxonomic tree.
Method \verb{nearest+} selects the traits belonging to the nearest taxa that have a taxonomic level equal or finer
than the target one. Method \verb{nearest-} does the opposite.
Method \verb{nearest+-} selects the traits belonging to the nearest taxa giving priority to taxa having
taxonomic level equal or finer than the target one. Method \verb{nearest-+} does the opposite.
}
\examples{
data(macro_ex)

data_bio <- as_biomonitor(macro_ex)
data_agr <- aggregate_taxa(data_bio)
data_ts <- assign_traits(data_agr)

# select only the nearest traits
data_ts_sub <- manage_traits(data_ts, method = "nearest+-")
}
