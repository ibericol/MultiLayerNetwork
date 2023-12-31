% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Networks Functions Final.R
\name{replace_with_tratified_mean}
\alias{replace_with_tratified_mean}
\title{Fills NA values with the mean by groups.}
\usage{
replace_with_tratified_mean(x, z)
}
\arguments{
\item{x}{a numeric vector with some NA.}

\item{z}{a numeric vector with group information. No NA alloewd.}
}
\value{
A numeric vector without NA.
}
\description{
`replace_with_mean_by_group()` Fills NA values with the mean by groups (specified by z). If
the group is missing, the NA value is replaced with the overall mean.
}
\examples{
replace_with_tratified_mean(x,z)
}
