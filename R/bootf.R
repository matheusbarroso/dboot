#'Garma statistic for tsboot
#'
#'\code{\link{bootf}} returns the estimated parameters of a GARMA
#'model and is used as an input (the \code{statistic}
#'parameter) for the \code{\link[boot]{tsboot}} function in
#' the \pkg{boot} package.
#'
#'The \emph{GARMA} model is estimated with a modified
#'version of the \code{\link[gamlss.util]{garmaFit}} function
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@examples
#'data(example_HHJ)
#'familia <- "PO"; ordem <- c(0,1)
#'bootf(example_HHJ)
#'
#'@param data A univariate or multivariate time series.
#' It might be vector, matrix or data frame to be passed
#'  to the \code{garmaFit2} function.
#'
#'@return A vector with the estimated Garma model parameters.


bootf <-
function (data) {
fit2 <- garmaFit2(yt~x-1,data=data,order=ordem,family=familia,tail=0,control=list(iter.max=1000))
return(fit2$coef)}
