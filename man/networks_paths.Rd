% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Networks Functions Final.R
\name{networks_paths}
\alias{networks_paths}
\title{Network Paths}
\usage{
networks_paths(W, transform = "Inv")
}
\arguments{
\item{W}{the Lower Level Network calculated with TOM.}

\item{transform}{Either Inv for inverse 1/W or Log for -log10 transformation.}
}
\value{
a list with distance matrix, the min distance matrix, a Next
matrix with the closest neighbor, and a path matrix with the best paths for e
ach pair of nodes.
}
\description{
`networks_paths()` calculates a distance matrix using the specified
transformation, then calculates the minimum distance matrix and the optimum
paths for each pair of nodes.
}
\examples{
split_path(path)
}
