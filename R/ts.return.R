#'Internal \code{\link{tsboot2}} function.
#'
#'Function to return the values \code{\link{tsboot2}}
#'
#'\code{ts.return} is an internal function called by
#'\code{\link{tsboot2}}. It is a modified version of the
#'\code{\link[boot]{ts.return}} function of \pkg{boot}.
#'
#'@author Matheus de Vasconcellos Barroso, \email{matheus.vb@gmail.com}
#'
#'@param t0 \code{t0} is the result of
#' \code{statistic(tseries,...{})}, as in \pkg{boot}
#'
#'@param t The results of applying \code{statistic} to the
#'replicate time series.
#'
#'@param R A positive integer giving the number of
#' bootstrap replicates required.
#'
#'@param tseries The original time series.
#'
#'@param seed The used seed in the bootstrap replication.
#'
#'@param stat \code{stat =  statistic}, as defined in
#' \code{\link{tsboot2}}
#'
#'@param n.sim The length of the simulated series.
#'
#'@param l \code{l} is the fixed block length used in generating the replicate
#'time series
#'
#'@param ran.gen This is a function of three arguments. The first
#' argument is a time series, it is the result of selecting
#' \code{n.sim} observations from \code{tseries} by some scheme
#' and converting the result back into a time series of the
#' same form as \code{tseries} (with length \code{n.sim}).
#' The second argument to \code{ran.gen} is always the value
#' \code{n.sim}, and the third argument is \code{ran.args},
#' which is used to supply any other objects needed by
#' \code{ran.gen}.
#'
#'@param ran.args This will be supplied to \code{ran.gen}
#'each time it is called. If \code{ran.gen} needs any
#'extra arguments then they should be supplied as components
#'of \code{ran.args}. Multiple arguments may be passed by
#'making \code{ran.args} a list. If \code{ran.args} is
#'\code{NULL} then it should not be used within
#'\code{ran.gen} but note that \code{ran.gen} must still
#'have its third argument.
#'
#'@param call The \code{\link{tsboot2}} call.
#'
#'@param blocks Blocks to be returned.
#'
#'@return A list with its arguments.
#'
#'@references
#'
#'@note For bugs and further requests please refer to
#' \url{https://github.com/matheusbarroso/dboot}


ts.return <-
function (t0, t, R, tseries, seed, stat, n.sim, l, ran.gen, ran.args, call, blocks) {
    out <- list(t0 = t0, t = t, R = R, data = tseries, seed = seed,sim='fixed',
        statistic = stat,  n.sim = n.sim, call = call,blocks=blocks,l = l)

        if (!is.null(call$ran.gen))
            out <- c(out, list(ran.gen = ran.gen, ran.args = ran.args))

    class(out) <- "boot"
    out
}
