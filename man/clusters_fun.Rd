% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Networks Functions Final.R
\name{clusters_fun}
\alias{clusters_fun}
\title{Cluster analysis on the distance matrix of TOM.}
\usage{
clusters_fun(W, ClGroup, ClType, GroupIndx_df)
}
\arguments{
\item{W}{the Lower Level Network calculated with TOM.}

\item{ClGroup}{group of variables that will be clustered.
If not specified, it will use all variables.}

\item{ClType}{one of "HC" (Hierarchical Clustering), "PAM"
(K Medoids clustering) or "NOCL" (No Cluster).}

\item{GroupIndx_df}{a data frame with the groups information.}
}
\value{
A list with the new group index, the silhouette plot, and the best K.
}
\description{
`clusters_fun()` does a cluster analysis on the TOM distance matrix
D and chooses the best number of clusters K by maximizing the silhouette
values for each possible number of clusters.
}
\examples{
clusters_fun(W,ClGroup,ClType,GroupIndx_df)
}
