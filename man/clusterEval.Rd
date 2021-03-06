% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{clusterEval}
\alias{clusterEval}
\title{Clustering Evaluation Against Ground Truth}
\usage{
clusterEval(x, y, verb = TRUE, perf = FALSE)
}
\arguments{
\item{x}{vector of size \code{n} describing a clustering.}

\item{y}{vector of size \code{n} describing a second clustering that will be compared to \code{x}.}

\item{verb}{boolean value: if \code{TRUE}, the function prints a line for each cluster.}

\item{perf}{boolean value: if \code{TRUE}, the function stops when it detects a different number of clusters and return 0.}
}
\value{
integer: percentage of correctly matched nodes when comparing the two clustering \code{x} and \code{y}.
}
\description{
This function evaluates the quality of a clustering when a ground truth clustering is available. It computes the
percentage of mismatched nodes when matching the clusters.
}
\examples{
g = generateEdgeList(n=1000)
res = ML_RMCL(g)
clusterEval(g$clusters, res)
}
\keyword{checking}
\keyword{clustering}
\keyword{evaluation}
\keyword{ground}
\keyword{truth}
