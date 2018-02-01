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
#'@param n.sim The length of the simulated series.
#'
#'@param seed The used seed in the bootstrap replication.
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
#'@param allow.parallel Logical TRUE/FALSE indicating
#' whether parallel computation via the foreach package
#' should be used. The default value is TRUE. OBS:paralllel
#'  backend must be registered prior to calling
#'  \code{\link{tsboot2}}.
#'
#'@param seed Numeric, the seed to \code{set.seed()} for
#' replicable examples.
#'
#'@param packages If \code{allow.parallel = TRUE}.
#' A character vector with the lisf of
#'packages required by \code{statisitc}.
#'
#'@param export If \code{allow.parallel = TRUE}.
#' A character vector with the lisf of
#'objects (functions, etc...) required by \code{statisitc}.
#'
#'@param \dots Extra argumetns to \code{statistic} may be
#'supplied here. Beware of partial matching to the
#'arguments of \code{\link{tsboot2}} listed above.
#'
#'@examples \dontrun{
#'library(dboot)
#'library(gamlss)
#'library(doParallel)
#'no_cores <- if(detectCores()==1) 1 else detectCores() -1
#'registerDoParallel(no_cores)
#'bootf <- function (db,ord,fam) {
#'  fit2 <- garmaFit2(yt~x-1,data=db,order=ord,family=fam,tail=0,control=list(iter.max=1000))
#'  return(fit2$coef)}
#'ord <- c(1,1) ; fam="GA"
#'db <- example_LFL
#'tsboot2(db, statistic=bootf, R = 10, l = 100,ord=ord,fam=fam,export=c("garmaFit2"),package=c("gamlss"))
#'}
#'
#'@return An object of class "boot" with similar components
#'as \code{\link[boot]{tsboot}}
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
function (tseries, statistic, R=100, l = NULL, n.sim = NROW(tseries),
		      ran.gen = function(tser,n.sim, args) tser, ran.args = NULL,
		      allow.parallel=TRUE, seed=123,packages=NULL,export=NULL, ...) {

  R <- floor(R)
  if((!is.numeric(R)||(R <=0)))
    stop("'R' must be positive integer")
  
  if(!is.logical(allow.parallel))
    stop("allow.parallel must be logical (i.e.: TRUE/FALSE)")
	
  if(allow.parallel)
		if (foreach::getDoParRegistered()==FALSE)
			stop("parallel backend must be registered")
	
  if(!is.numeric(seed))
	  stop('seed must be numeric')

  if(!is.null(packages)&&!is.character(packages))
    stop('packages must be a character vector')

  if(!is.null(export)&&!is.character(export))
    stop('packages must be a character vector')

	statistic
  tscl <- class(tseries)
	call <- match.call()
  t0 <- statistic(tseries, ...)
  
  ts.orig <- if (!boot:::isMatrix(tseries))
    as.matrix(tseries)
      else tseries
  
  n <- nrow(ts.orig)
  
  if (missing(n.sim))
    n.sim <- n
  
  class(ts.orig) <- tscl
  
    if ((is.null(l) || (l <= 0) || (l > n)))
        stop("invalid value of 'l'")

    fn <-  {
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

	`%op%` <- if(allow.parallel==TRUE) `%dorng%` else `%do%`
	set.seed(seed) # or .options.RNG=123; mudar para inherit
	res <- foreach(i=seq_len(R),.packages=packages,.export=export)%op%fn(i)
    t <- matrix(, R, length(res[[1L]]$statistic))
    for (r in seq_len(R)) t[r, ] <- res[[r]]$statistic
    ts.return(t0 = t0, t = t, R = R, tseries = tseries, seed = seed,
        stat = statistic, n.sim = n.sim,
        l = l, ran.gen = ran.gen, ran.args = ran.args, call = call,
        blocks=res[[1L]]$blocks)
}
