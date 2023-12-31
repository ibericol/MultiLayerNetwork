% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Networks Functions Final.R
\name{get_path}
\alias{get_path}
\title{Get Path}
\usage{
get_path(Next, start, end)
}
\arguments{
\item{Next}{a Matrix with the next closest neighbor between two nodes.}

\item{start}{the starting node.}

\item{end}{the ending node.}
}
\value{
a string with the best path from the starting node to the ending node.
}
\description{
`get_path()` gets the best path between two nodes.
}
\examples{
floyd_warshall(Distance)
}
