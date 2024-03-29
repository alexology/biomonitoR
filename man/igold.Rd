% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/igold.R
\name{igold}
\alias{igold}
\title{igold}
\usage{
igold(x, traceB = FALSE)
}
\arguments{
\item{x}{Results of \code{aggregate_taxa()}.}

\item{traceB}{If set to \code{TRUE} a list as specified below will be returned.}
}
\value{
If \code{traceB} is set to \code{TRUE} a list with the following elements will be returned:
\itemize{
\item \code{results} Results of the \code{igold} index.
\item \code{taxa_df} The \code{data.frame} used for the calculation containing the abundance of the GOLD taxa.
}
}
\description{
This function calculates the 1 - GOLD metric, where GOLD stands for Gastropoda, Oligochaeta and Diptera. This metric should decrease with increasing organic pollution (Pinto et al., 2004).
}
\details{
The metric 1 - GOLD is calculated as 1 minus the relative abundance of Gastropoda, Oligochaeta and Diptera. If a custom database is provided (see \link{aggregate_taxa}) please be sure that Gastropoda are at Class, Oligochaeta at Sublclass and Diptera at Order level, otherwise the gold calculation will be meaningless.
If this is the case please see \code{get_taxa_abundance()}.
}
\examples{
data(macro_ex)
data_bio <- as_biomonitor(macro_ex)
data_agr <- aggregate_taxa(data_bio)
igold(data_agr)
}
\references{
Pinto, P., Rosado, J., Morais, M., & Antunes, I. (2004). Assessment methodology for southern siliceous basins in Portugal. In Integrated Assessment of Running Waters in Europe (pp. 191-214). Springer, Dordrecht.
}
\seealso{
\code{\link{aggregatoR}}
}
\keyword{ept}
