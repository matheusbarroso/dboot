#'LFL-MBB optimal block lenght selection algorithm
#'
#'This functions performs the data-based algorithm of Lahiri,
#' Furukawa and Lee (2005), henceforth \strong{LFL}, for
#'  the selection of optimal block length sizes in the case
#'   of block bootstrap of Kunsch (1989).
#'
#'This functions implements the iterative version of the
#'Lahiri, Furukawa and Lee (2005) algorithm. Here some
#'modifications are implemented in the fashion of Barroso
#'(2017). Namely, a vectorized algorithm is implemented
#'where the user might supply which parameter to optimize
#'over or use a default value. The default value is obtained
#' by minimizin the mean MSE vector (if a vector or
#' parameters is returned by statistic). For the MBB step
#' the \code{\link{tsboot2}},a modification of the
#' \code{\link[boot]{tsboot}} function from the
#' \pkg{boot} package, is used.
#'
#'@author Matheus de Vasconcellos Barroso
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
#'
#'LFL(db,bootf,ord=ord,fam=fam,export=c("garmaFit2"),package=c("gamlss"))}

#'@param data A univariate or multivariate time series.
#'It might be vector, matrix or data frame to be passed to
#'statistic.
#'
#'@param statistic A function which when applied to data
#' returns a vector containing the statistic of interest.
#' Each time statistic is called it is passed a time series
#'  of length n which is of the same class as the original
#'  tseries. Any other arguments which statistic takes must
#'   remain constant for each bootstrap replicate and
#'   should be supplied through the . . . argument to
#'   tsboot.
#'
#'@param R A positive integer giving the number of
#' bootstrap replicates required. The default value is 100.
#'
#'@param nsteps A positive integer with the number of steps
#'(iterations) in the HHJ algorithm. The default value is 5.
#'
#'@param l.init A positive integer smaller then the number
#' of observations (n.obs.) in \code{data} (i.e. the length
#'  of the time series), indicating the size of the subset
#'  in the \bold{HHJ} algorithm. The default value is
#'  'default', in which \eqn{m.init = round(c1*n.obs.^(1/(r+4)))}
#'
#'@param type.optm 0 for mean of the parameters vector
#' or a positive integer giving the index of the desired
#' parameter to optimize in the same order as provided by
#' statistic. For more details see the section
#' \strong{"Details"} bellow.#'
#'
#'@param type.est A character describing the type of
#'estimation being undertaken. Accepted values are:
#'"bias.variance" and "distribution.quantile".
#' The default value is 'bias.variance'.
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
#'  backend must be registered prior to calling \strong{HHJ}.
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
#'
#'@return A dataframe, the first column for the iteration
#'step the second and third for the estimated optimal block length.
#'
#'@references Kunsch, Hans R. 1989. The Jakknife and the
#'Bootstrap For General Stationary Observations.
#'The Annals of Statistics. 1989, Vol. 17, 3, pp. 1217-1241.
#'
#'Lahiri, Soumendra, Furukawa, Kyoji and Lee, Yd. 2007.
#'A nonparametric plug-in rule for selecting optimal block lengths
#'for block bootstrap methods. Statistical Methodology. July, 2007,
#'Vol. 4, 3, pp. 292-321.
#'
#'Barroso, Matheus de V. 2018.  BOOTSTRAP METHODS FOR
#'GENERALIZED AUTOREGRESSIVE MOVING AVERAGE MODELS
#'
#'@note This function is not working properly. Wait for an
#'update.
#'For bugs and further requests please refer
#' to \url{https://github.com/matheusbarroso/dboot}



LFL <-
function(data,statistic,R=100L,nsteps=5L,
		 l.init=NULL, type.optm=0,
		 type.est="bias.variance",
		 ran.gen = function(tser, n.sim, args) tser,
		 ran.args=NULL,allow.parallel = TRUE,
		 seed = 123, packages = NULL, export = NULL,...) {

R <- floor(R)
if((!is.numeric(R)||(R <=0)))
	stop("'R' must be positive integer")

nsteps <- floor(nsteps)
if (!(is.numeric(nsteps)&(nsteps > 0 )))
	stop("nsteps must be a positive integer ")

if(!is.numeric(type.optm)||(type.optm <0))
	stop("type.otpm must be a positive integer; accepted values are:
         0 for mean of the parameters vector or the index of the desired parameter
         to optimize in the same order as provided by statistic")

if (type.est%in% c("bias.variance","distribution.quantile") ) {
	switch(type.est,
	"bias.variance"={
		r <- 1
		a <- 0
	     	c2 <- 1
			}
	 ,"distribution.quantile"={
		r <- 2
		a <- 1/2
	    	c2 <- 0.1
			          }
		  )
} else stop("unrecognised value of 'type.est',accepted values are:
			'bias.variance' and 'distribution.quantile'")

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

n <- if(is.vector(data)) length(data) else
		if(is.matrix(data)||is.data.frame(data)) nrow(data)

c1 <- 1

if(is.null(l.init)) {
	l.init <- c1*n^(1/(r+4))
	l.init <- ceiling(l.init)
					}

opt.l<- list(iteration=seq.int(from=0,nsteps),l=list(),l.adj=list())
opt.l$l[[1]] <- l.init
opt.l$l.adj[[1]] <- l.init

for (iteration in seq_len(nsteps)) {

	##step 1
	l <- if(type.optm==0) ceiling(mean(opt.l$l.adj[[iteration]])) else {
			if(is.numeric(type.optm)&(type.optm > 0))
				opt.l$l.adj[[iteration]][type.optm] else
					stop("'type.optm' must be a positive integer,smaller than or
						 equal to the number of parameters in the model.")
																		}
	if(l < 1) {
		l <- 1L
		print("The ajusted value of l was smaller than one and was replaced by one/n")
			  }

	estimate.l <- tsboot2(data, statistic = statistic, R = R, l = l,ran.gen = ran.gen,
						  ran.args = ran.args,allow.parallel = allow.parallel,
						  seed = seed, packages = packages, export = export,...)

	estimate.2l <- tsboot2(data, statistic = statistic, R = R, l = 2*l, ran.gen = ran.gen,
						  ran.args = ran.args,allow.parallel = allow.parallel,
						  seed = seed, packages = packages, export = export,...)

	print(paste("Finished the boot step",iteration))
	##step 2
	endpt <- n - l + 1  ## endpoint
	m <- c2*n^(1/3)*l^(2/3)  ## m= number of bootstrap blocks to be deleted
	m <- ceiling (m)
	if(m >= endpt) stop("invalid value of m: m < n-3+1")
	M <- endpt -m +1
	K <- 1:R
	group <- M.collection(m=m,M=M,R=R,blocks=estimate.l$blocks) ## parallelizar

	block.jab.point.value <- lapply(group, function(k) {

		if (is.null(k)) block.jab.point.value <- NULL else {

			x <- estimate.l$blocks$starts[k,]
			if(is.null(dim(x))) {

				blocks <- sapply(seq_len(dim(estimate.l$blocks$starts)[1]),
				function(j) {all(x%in%estimate.l$blocks$starts[j,])})

				block.jab.point.value <-  estimate.l$t[blocks,]
								} else {

					blocks <- apply(x,1,function (i) {sapply(seq_len(dim(estimate.l$blocks$starts)[1]),
					function(j) {all(i%in%estimate.l$blocks$starts[j,])})})

					boot.blocks <- apply(blocks,2,function(j) estimate.l$t[j,])

					block.jab.point.value <- apply(boot.blocks,1,mean)
										}
															}
														})

	block.jab.point.value <-matrix(unlist(block.jab.point.value),ncol=dim(estimate.l$t)[2],byrow=T)

	pseudo <- (matrix(estimate.l$t0,ncol=dim(block.jab.point.value)[2],
	nrow=dim(block.jab.point.value)[1],byrow=T)*endpt-(endpt-m)*block.jab.point.value)/m

	var.jab <- apply(apply(pseudo,1,function(i)(i-estimate.l$t0)^2 ),1,sum)*m/((endpt-m)*M)

	##step 3
	C1.hat <- n*l^(-r)*var.jab
	C2.hat <- 2*l*(estimate.l$t0-estimate.2l$t0)
	l.optm <- (2*C2.hat^2/(r*C1.hat))^(1/(r+2))*n^(1/(r+2))
	opt.l$l[[iteration+1]] <- l.optm
	opt.l$l.adj[[iteration+1]] <- l.optm^((r+2)/(r+4))

	print(paste("Step:",iteration, "completed, more", nsteps-iteration,"to go\n"))
	print(paste("Initial l is:",l.init,"Estimated length in step",iteration,": original = ",
		opt.l$l[[iteration+1]],", adjusted = ",opt.l$l.adj[[iteration+1]],"\n"))

}
return(opt.l)

}



