% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{summarise_num}
\alias{summarise_num}
\title{Summarise numeric vector}
\usage{
summarise_num(x)
}
\arguments{
\item{x}{a numeric vector.}
}
\value{
A numeric vector with five named elements: \code{mean}, \code{sd}, \code{median},
\code{p25}, and \code{p75}, corresponding to the mean, standard deviation, median,
and 25-/75-percentiles.
}
\description{
Used internally, to summarise numeric vectors.
}
\keyword{internal}
