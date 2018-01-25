#'Garma statistic for tsboot
#'
#'This function returns the estimated parameters of a GARMA
#'model and is used as an input (the \code{statistic}
#'parameter) for the \code{\link[boot]{tsboot}} function in the \pkg{boot} package.
#'
#'@author Matheus de Vasconcellos Barroso
#'@example
#'data(example_HHJ)
#'familia <- "PO"; ordem <- c(0,1)
#'bootf(example_HHJ)
#'@param data A univariate or multivariate time series.
#' It might be vector, matrix or data frame to be passed
#'  to the \code{garmaFit2} function.?

bootf <-
function (data) {
fit2 <- garmaFit2(yt~x-1,data=data,order=ordem,family=familia,tail=0,control=list(iter.max=1000))
return(fit2$coef)}
