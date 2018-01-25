#'Internal \code{\link{LFL}} function.
#'
#'Generate R bootstrap replicates of a statistic applied
#'to a time series. The replicate time series can be
#'generated using fixed lengths replicates.
#'
#'This is a modified version of the
#'\code{\link[boot]{tsboot}} of the \pkg{boot} package.
#'Here, only Moving Block Bootstrap is implemented,
#'the function is modified to return the resampled blocks
#' and some unecessary arguments were deleted. Additionaly,
#' the \code{parallel} option is through
#' the \pkg{doParallel} package.
#'
#'@author Matheus de Vasconcellos Barroso, \email{matheus.vb@gmail.com}
#'
#'@param tseries A univariate or multivariate time series.
#'
#'@param statistic A function which when applied to
#'\code{tseries} returns a vector containing the
#'statistic(s) of interest. Each time \code{statistic} is
#'called it is passed a time series of length \code{n.sim}
#'which is of the same class as the original \code{tseries}.
#'Any other arguments which \code{statistic}  takes must
#'remain constant for each bootstrap replicate and should
#'be supplied through the \dots argument to
#'\code{\link{tsboot2}}.
#'
#'@param R A positive integer giving the number of
#' bootstrap replicates required.
#'
#'@param l \code{l} is the fixed block length used in generating the replicate
#'time series
#'
#'@param sim A character with the 'fixed' tag. (REMOVE the par.)
#'
#'@param endcorr \code{FALSE} (REMOVE the par.)
#'
#'@param n.sim The length of the simulated series. (REMOVE the par.)

#'
#'@param seed The used seed in the bootstrap replication.
#'
#'@param orig.t (REMOVE the par.)
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
#'@param norm Not used (REMOVE the par.)
#'
#'@param parallel Not used (REMOVE the par.)
#'
#'@param ncpus Not used (REMOVE the par.)
#'
#'@param cl Not used (REMOVE the par.)
#'
#'@allow.parallel Logical TRUE/FALSE indicating
#' whether parallel computation via the foreach package
#' should be used. The default value is TRUE. OBS:paralllel
#'  backend must be registered prior to calling
#'  \code{\link{tsboot2}}.
#'
#'@param \dots Blocks to be returned.
#'

#'
#'@return A list with the two elements:
#' \item{starts }{The initial position of the resampled blocks}
#' \item{lengths }{The length of the blocks in \code{starts}}
#'
#'@references
#'
#'Angelo Canty and Brian Ripley (2017). boot: Bootstrap R
#'(S-Plus) Functions. R package version 1.3-19.
#'
#'Davison, A. C. & Hinkley, D. V. (1997) Bootstrap
#'Methods and Their Applications. Cambridge University
#'Press, Cambridge. ISBN 0-521-57391-2
#'
#'@note For bugs and further requests please refer to
#' \url{https://github.com/matheusbarroso/dboot}



tsboot2 <-
function (tseries, statistic, R, l = NULL, sim = "fixed", endcorr = TRUE,
    n.sim = NROW(tseries), orig.t = TRUE, ran.gen = function(tser,
        n.sim, args) tser, ran.args = NULL, norm = TRUE, ...,
    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus",
        1L), cl = NULL) {
    if (missing(parallel))
        parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
    if (parallel != "no" && ncpus > 1L) {
        if (parallel == "multicore")
            have_mc <- .Platform$OS.type != "windows"
        else if (parallel == "snow")
            have_snow <- TRUE
        if (!have_mc && !have_snow)
            ncpus <- 1L
        loadNamespace("parallel")
    }
    statistic
    tscl <- class(tseries)
    R <- floor(R)
    if (R <= 0)
        stop("'R' must be positive")
    call <- match.call()
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    t0 <- if (orig.t)
        statistic(tseries, ...)
    else NULL
    ts.orig <- if (!boot:::isMatrix(tseries))
        as.matrix(tseries)
    else tseries
    n <- nrow(ts.orig)
    if (missing(n.sim))
        n.sim <- n
    class(ts.orig) <- tscl
    if ((is.null(l) || (l <= 0) || (l > n)))
        stop("invalid value of 'l'")
i.a.2 <- list()
    fn <- if (sim %in% c("fixed")) {
        i.a <- ts.array(n, n.sim, R, l)
        ran.gen
        ran.args
        function(r) {
            ends <-  cbind(i.a$starts[r, ], i.a$lengths)
            inds <- apply(ends, 1L, boot:::make.ends, n)
            inds <- if (is.list(inds))
                matrix(unlist(inds)[1L:n.sim], n.sim, 1L)
            else matrix(inds, n.sim, 1L)
            list(statistic=statistic(ran.gen(ts.orig[inds, ], n.sim, ran.args), ...),blocks=i.a)

        }
    }
    else stop("unrecognized value of 'sim'")
    res <- if (ncpus > 1L && (have_mc || have_snow)) {
	stop("parallel computation not implemented, yet") ##não está implementado... tem que alterar o código.
        if (have_mc) {
            parallel::mclapply(seq_len(R), fn, mc.cores = ncpus)
        }
        else if (have_snow) {
            list(...)
            if (is.null(cl)) {
                cl <- parallel::makePSOCKcluster(rep("localhost",
                  ncpus))
                if (RNGkind()[1L] == "L'Ecuyer-CMRG")
                  parallel::clusterSetRNGStream(cl)
                res <- parallel::parLapply(cl, seq_len(R), fn)
                parallel::stopCluster(cl)
                res
            }
            else parallel::parLapply(cl, seq_len(R), fn)
        }
    }
    else lapply(seq_len(R), fn)
    t <- matrix(, R, length(res[[1L]]$statistic))
    for (r in seq_len(R)) t[r, ] <- res[[r]]$statistic
    ts.return(t0 = t0, t = t, R = R, tseries = tseries, seed = seed,
        stat = statistic, sim = sim, endcorr = endcorr, n.sim = n.sim,
        l = l, ran.gen = ran.gen, ran.args = ran.args, call = call,
        norm = norm,blocks=res[[1L]]$blocks)
}
