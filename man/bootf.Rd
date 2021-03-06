% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootf.R
\name{bootf}
\alias{bootf}
\title{Garma statistic for tsboot}
\usage{
bootf(data, ord, fam)
}
\arguments{
\item{data}{A univariate or multivariate time series.
It might be vector, matrix or data frame to be passed
 to the \code{garmaFit2} function.}
}
\value{
A vector with the estimated Garma model parameters.
}
\description{
\code{\link{bootf}} returns the estimated parameters of a GARMA
model and is used as an input (the \code{statistic}
parameter) for the \code{\link[boot]{tsboot}} function in
the \pkg{boot} package.
}
\details{
The \emph{GARMA} model is estimated with a modified
version of the \code{\link[gamlss.util]{garmaFit}} function
}
\examples{
example_HHJ
fam <- "PO"; ord <- c(0,1)
bootf(example_HHJ,ord,fam)

}
\author{
Matheus de Vasconcellos Barroso
}
