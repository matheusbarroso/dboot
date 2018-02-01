#'HHJ-MBB optimal block lenght selection algorithm
#'
#'This functions performs the data-based algorithm of Hall,
#' Horowitz and Jing (1995), henceforth \strong{HHJ}, for
#'  the selection of optimal block length sizes in the case
#'   of block bootstrap of Kunsch (1989).
#'
#'This functions implements the iterative version of the
#'Hall, Horowitz and Jing (1995) algorithm. Here some
#'modifications are implemented in the fashion of Barroso
#'(2017). Namely, a vectorized algorithm is implemented
#'where the user might supply which parameter to optimize
#'over or use a default value. The default value is obtained
#' by minimizin the mean MSE vector (if a vector or
#' parameters is returned by statistic). For the MBB step
#' the \code{\link[boot]{tsboot}} function from the
#' \pkg{boot} package is used. Additionaly, with the
#'  argument \code{type.sub.blocks} the user might use a
#'  sampling sche on the subsets of \code{m.init}. There is
#'   one sample from the quantile 1/3, another for the 2/3
#'   and the last one for the 1/3-2/3.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#'@examples \dontrun{
#'library(gamlss)
#'library(dboot)
#'db <- example_HHJ
#'fam <- "PO"
#'ord <- c(0,1)
#'
#'library(doParallel)
#'no_cores <- if(detectCores()==1) 1 else detectCores() -1
#'registerDoParallel(no_cores)
#'
#'HHJ(db,bootf,R = 10,nsteps=1,ord=ord,fam=fam,
#'export=c("garmaFit2"), packages=c("gamlss"))
#'
#'HHJ(db,bootf,R = 10, errorhandling="pass",ord=ord,seed=124
#'fam=fam,export=c("garmaFit2"),  packages=c("gamlss"),n.try=3)
#'
#'HHJ(db,bootf,R = 10,nsteps=2,ord=ord,fam=fam,
#'export=c("garmaFit2"),  packages=c("gamlss"),n.try=2)
#'
#'HHJ(db,bootf,R = 10,nsteps=1,ord=ord,fam=fam,
#'export=c("garmaFit2"), packages=c("gamlss"),type.sub.blocks="complete")
#'
#'HHJ(db,bootf,R = 10,nsteps=2,ord=ord,fam=fam,
#'export=c("garmaFit2"),  packages=c("gamlss"),type.optm=2)
#'
#'HHJ(db,bootf,R = 10,nsteps=3,ord=ord,fam=fam,
#'export=c("garmaFit2"), packages=c("gamlss"),m.init=6)
#'
#'HHJ(db,bootf,R = 10,nsteps=1,ord=ord,fam=fam,
#'export=c("garmaFit2"), packages=c("gamlss"),m.init=6,allow.parallel=FALSE)
#'
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
#'@param type.est A character describing the type of
#'estimation being undertaken. Accepted values are:
#''bias.variance', 'one.sided.distribution' and
#' 'two.sided.distribution. The default value is
#'  'bias.variance'.
#'
#'@param type.optm 0 for mean of the parameters vector
#' or a positive integer giving the index of the desired
#' parameter to optimize in the same order as provided by
#' statistic. For more details see the section
#' \strong{"Details"} bellow.
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
#'@param errorhandling A character, accepted values are 'try' and 'pass'.
#'
#'@param n.try A positive integer indicating the maximum
#' number of tries if errorhandling='try'. The default
#' value is 3.
#'
#'@param m.init A positive integer smaller then the number
#' of observations (n.obs.) in \code{data} (i.e. the length
#'  of the time series), indicating the size of the subset
#'  in the \bold{HHJ} algorithm. The default value is
#'  'default', in which \code{m.init=round(n.obs.^(1/3))}
#'
#'@param type.sub.blocks A character, accepted values are: 'fast'
#' and 'complete', the former implements a sampling scheme
#'  in the subsets of \code{m (or m.init)}, where only
#'  three subsets lengths are evaluated. In the latter all
#'  subsets lenghts are evaluated. For a large dataset the
#'  option 'complete' might take too long. For more details
#'   see the section \strong{"Details"} bellow.
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
#'arguments of \code{\link{tsboot2}}.
#'
#'@return A dataframe, the first column for the iteration
#'step and the second for the estimated optimal block length.
#'
#'@references Kunsch, Hans R. 1989. The Jakknife and the
#'Bootstrap For General Stationary Observations.
#'The Annals of Statistics. 1989, Vol. 17, 3, pp. 1217-1241.
#'
#'Hall, P., Horowitz, J. L e Jing, B-Y. 1995. On blocking
#'rules for the bootstrap with dependent data. Biometrika. 82, 1995, pp. 561-574.
#'
#'Barroso, Matheus de V. 2018.  BOOTSTRAP METHODS FOR
#'GENERALIZED AUTOREGRESSIVE MOVING AVERAGE MODELS
#'
#'@note For bugs and further requests please refer
#' to \url{https://github.com/matheusbarroso/dboot}





HHJ <-
function(data, statistic, R = 100L, nsteps = 5L,
		 type.est = "bias.variance", type.optm = 0,
		 ran.gen = function(tser, n.sim, args) tser,
		 ran.args = NULL,allow.parallel = TRUE,
		 errorhandling = "try",n.try = 5L,
		 m.init = "default",type.sub.blocks = "fast",
		 seed = 123, packages = NULL,export=NULL,...) {

R <- floor(R)
if((!is.numeric(R)||(R <=0)))
	stop("'R' must be positive integer")

nsteps <- floor(nsteps)
if((!is.numeric(nsteps)||(nsteps <=0)))
	stop("nsteps must be a positive integer")

if (type.est%in% c("bias.variance","one.sided.distribution","two.sided.distribution") ) {
	switch(type.est,"bias.variance"={k <- 1}
				   ,"one.sided.distribution"={k <- 2}
				   ,"two.sided.distribution"={k <- 3}
		  )
}else stop("unrecognised value of 'type.est',accepted values are: 'bias.variance',
		   'one.sided.distribution' and 'two.sided.distribution'")

if(!is.numeric(type.optm)||(type.optm <0))
	stop("type.otpm must be a positive integer; accepted values are:
	0 for mean of the parameters vector or the index of the desired parameter
	to optimize in the same order as provided by statistic")

if(!is.logical(allow.parallel))
	stop("allow.parallel must be logical (i.e.: TRUE/FALSE)")

if(allow.parallel)
	if (foreach::getDoParRegistered()==FALSE)
		stop("parallel backend must be registered")

if (!errorhandling%in% c("try","pass") )
	stop("unrecognised value of 'errorhandling',accepted values are: 'try' and 'pass'")

if((errorhandling=="try")&&(!is.numeric(n.try)||(n.try <=0)))
	stop("n.try must be a positive integer")

n <- if(is.vector(data)) length(data) else 
	if(is.matrix(data)||is.data.frame(data)) nrow(data)

if(m.init!="default")
	if(!(is.numeric(m.init)&&(m.init < n)&&(m.init > 0)))
		stop("m.init must be a positive integer smaller than n. \nTip: you could also specify m.init='default' and in this case m.init=round(n^(1/3))")

if ( !type.sub.blocks%in% c("fast","complete") )
	stop("unrecognised value of 'type.sub.blocks',accepted values are: 'fast' and 'complete'")

if(!is.numeric(seed))
	stop('seed must be numeric')

if(!is.null(packages)&&!is.character(packages))
	stop('packages must be a character vector')

if(!is.null(export)&&!is.character(export))
	stop('packages must be a character vector')	

	
opt.l <- data.frame(iter=seq.int(from=0,to=nsteps),l.optm=c(if(m.init=="default")round(n^(1/3)) else 
	round(m.init) ,rep(0,nsteps)))


set.seed(seed)
t0 <- statistic(data,...)

for(iteration in seq_len(nsteps)) {

	m <- opt.l$l.optm[iteration] # m = subsample size
	m.subsamples <- lapply(0:(n-m), function(j) {1:m+j}) ## all blocks of length m
	block.l.m.sub <- seq_len(m) 	## possible block lengths in m.subsamples
	if(type.sub.blocks=="fast") {
		inf <- seq_len(quantile(block.l.m.sub,1/3,type=1))
		sup <- seq.int(quantile(block.l.m.sub,2/3,type=1),m)
		dif <- setdiff(block.l.m.sub ,c(inf,sup))
		block.l.m.sub <- c(sample(inf,1),if(!identical(dif,integer(0))) sample(dif,1),sample(sup,1))
		block.l.m.sub <- unique(block.l.m.sub)
								}
	cat(paste("Chosen block lengths:",block.l.m.sub,"\n"))
    cat(paste("Wait while the MBB resamples for the subsample are being computed;",type.sub.blocks,"computation enabled\n"))

	h <- function(w) if( any( grepl( "<anonymous>:...", w) ) ) invokeRestart( "muffleWarning" )

	
	fit.mbb <- { 
		lapply(block.l.m.sub, function(l) {
			lapply(m.subsamples, function(j){
			switch(errorhandling
				,"try"= {
					MBB <- NULL
					attempt <- 1
					while( is.null(MBB) && attempt <= n.try ) {
					attempt <- attempt + 1
					try(MBB <- withCallingHandlers( tsboot2(subset(data,seq_len(n)%in%j), statistic = statistic, 
										R = R, l = l,ran.gen = ran.gen,
										ran.args = ran.args,allow.parallel = allow.parallel,
										seed = seed + attempt -1, packages = packages, export = export,...), 
										warning = h ) ) 
															  }
														   
														   
						}
				,"pass"={MBB <- withCallingHandlers( tsboot2(subset(data,seq_len(n)%in%j), statistic = statistic,
										R = R, l = l,ran.gen = ran.gen,
										ran.args = ran.args,allow.parallel = allow.parallel,
										seed = seed, packages = packages, export = export,...), 
										warning = h )})

			 list(MBB=MBB,seed=seed,subsamble=j,
			 index=seq_len(length(m.subsamples))[unlist(lapply(m.subsamples, function(subsamp) all(subsamp==j)))])
			
												})
											})
			
				}

	check01 <- lapply(seq_len(length(block.l.m.sub)), function(k)setdiff(seq_len(length(m.subsamples)),sapply(fit.mbb[[k]],function(j) j$index)) )
	if(!identical(unlist(check01),integer(0))) warning(sapply(seq_len(length(block.l.m.sub)), function(j) paste("length",block.l.m.sub[j],", iterations ",
	paste(as.character(if(identical(check01[[j]],integer(0)))0 else check01[[j]]),collapse=",")," have failed in step",iteration,";\n " ,sep="")))

	check <- lapply(seq_len(length(block.l.m.sub)), function(k) lapply(fit.mbb[[k]],function(j) if(is.null(j$MBB)) j$index))
	if(!identical(unlist(check),NULL)) warning(sapply(seq_len(length(block.l.m.sub)), function(j) paste("length",j,", iterations ",
	paste(sapply(check,unlist)[[j]],collapse=",")," have failed in step",iteration,";\n" ,sep="")))

	index1 <- block.l.m.sub[sapply(fit.mbb,length)!=0] ## for a given block length there might be a failure of the MBB
	index2 <- block.l.m.sub[!sapply(fit.mbb,is.null)]
	index1 <- seq_len(length(block.l.m.sub))[block.l.m.sub%in%intersect(index1,index2)]
	
	
	sqe <- lapply(index1, function(k) {
	sqe <- sapply(seq_len(length(fit.mbb[[k]]))[setdiff(seq_len(length(fit.mbb[[k]])),unlist(check[[k]]))], function(j) {
	tn <- apply(fit.mbb[[k]][[j]]$MBB$t,2,mean)
	(t0-tn)^2})})
	
	(msqe <- lapply(sqe,function(j) {apply(j,1,mean)}))

	indice.best <- sapply(seq_len(length(msqe)), function(j) {mean(msqe[[j]][if(type.optm==0) TRUE else type.optm ])})
	l.m <- index1[which.min(indice.best)]

	opt.l$l.optm[iteration+1] <- round((n/length(block.l.m.sub))^(1/3)*l.m)

	cat(paste("Step:",iteration, "completed, more", nsteps-iteration,"to go\n"))
}
return(opt.l)
}