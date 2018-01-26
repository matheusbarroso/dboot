#'Internal \code{\link{tsboot2}} function.
#'
#'Computes the resampled blocks initial position and length.
#'
#'\code{ts.array} is an internal function called by
#'\code{\link{tsboot2}}. It is a modified version of the
#'\code{\link[boot]{ts.array}} function of \pkg{boot}.
#'
#'@author Matheus de Vasconcellos Barroso, \email{matheus.vb@gmail.com}
#'
#'@param n \code{n} the length of the bootstraped time series.
#'\eqn{n = nrow(tseries)}.
#'
#'@param n.sim The length of the simulated series.
#'
#'@param R A positive integer giving the number of
#' bootstrap replicates required.
#'
#'@param l \code{l} is the fixed block length used in generating the replicate
#'time series
#'
#'@return A list with the two elements:
#' \item{starts }{The initial position of the resampled blocks}
#' \item{lengths }{The length of the blocks in \code{starts}}
#'
#'@references
#'
#'Lahiri, S. N. 2002. ON THE JACKKNIFE-AFTER-BOOTSTRAP METHOD
#' FOR DEPENDENT DATA AND ITS CONSISTENCY PROPERTIES.
#' Econometric Theory. 18, 2002, pp. 79-98.
#'
#'Lahiri, Soumendra, Furukawa, Kyoji and Lee, Yd. 2007.
#' A nonparametric plug-in rule for selecting optimal block
#' lengths for block bootstrap methods. Statistical
#' Methodology. July, 2007, Vol. 4, 3, pp. 292-321.
#'
#'Barroso, Matheus de V. 2018.  BOOTSTRAP METHODS FOR
#'GENERALIZED AUTOREGRESSIVE MOVING AVERAGE MODELS
#'
#'@note For bugs and further requests please refer to
#' \url{https://github.com/matheusbarroso/dboot}



ts.array <-
function (n, n.sim, R, l) {
    endpt <- n - l + 1 ##total number of observer blocks of length l
    cont <- TRUE
    {
        nn <- ceiling(n.sim/l)
        lens <- c(rep(l, nn - 1), 1 + (n.sim - 1)%%l)
        st <- matrix(sample.int(endpt, nn * R, replace = TRUE),
            R)
    }
    list(starts = st, lengths = lens)
}
