HHJ <-
function(data,statistic,R=100,nsteps=5L,type.est="bias.variance",type.optm=0,
		 allow.parallel=TRUE,errorhandling="try",n.try=3,m.init="default",
		 type.sub.blocks="fast",seed=123,export,...) { 

if((!is.numeric(nsteps)||(nsteps <=0)))
	stop("nsteps must be a positive integer")
	 
if(!is.logical(allow.parallel)) 
	stop("allow.parallel must be logical (i.e.: TRUE/FALSE)")

if(allow.parallel) 
	if (foreach::getDoParRegistered()==FALSE) 
		stop("paralllel backend must be registered")

if (! errorhandling%in% c("try","pass") ) 
	stop("unrecognised value of 'errorhandling',accepted values are: 'try' and 'pass'")

if((errorhandling=="try")&&(!is.numeric(n.try)||(n.try <=0)))   
	stop("n.try must be a positive integer")

if(!is.numeric(type.optm)||(type.optm <0))   
	stop("type.otpm must be a positive integer; accepted values are: 
	0 for mean of the parameters vector or the index of the desired parameter to optimize in the same order as provided by statistic")

n <- if(is.vector(data)) length(data) else if(is.matrix(data)||is.data.frame(data)) nrow(data)

if(m.init!="default") 
	if(!(is.numeric(m.init)&&(m.init < n)&&(m.init > 0))) 
		stop("m.init must be a positive integer smaller than n. \nTip: you could also specify m.init='default' and in this case m.init=round(n^(1/3))")

if(!is.numeric(seed)) 
	stop('seed must be numeric')

if (type.est%in% c("bias.variance","one.sided.distribution","two.sided.distribution") ) {
	switch(type.est,"bias.variance"={k <- 1}
				   ,"one.sided.distribution"={k <- 2}
				   ,"two.sided.distribution"={k <- 3}							
		  )
}else stop("unrecognised value of 'type.est',accepted values are: 'bias.variance', 'one.sided.distribution' and 'two.sided.distribution'")

if ( !type.sub.blocks%in% c("fast","complete") ) 
	stop("unrecognised value of 'type.sub.blocks',accepted values are: 'fast' and 'complete'")

opt.l <- data.frame(iter=seq.int(from=0,to=nsteps),l.optm=c(if(m.init=="default")round(n^(1/3)) else round(m.init) ,rep(0,nsteps)))

`%op%` <- if(allow.parallel==TRUE) `%dopar%` else `%do%`

set.seed(seed)

t0 <- statistic(data)

for(iteration in seq_len(nsteps)) { 

	m <- opt.l$l.optm[iteration] # m = subsample size
	m.subsamples <- lapply(0:(n-m), function(j) {1:m+j}) ## all blocks of length m
	block.l.m.sub <- seq_len(m) 	## possible block lengths in m.subsamples
	if(type.sub.blocks=="fast") {
		inf <- seq_len(quantile(block.l.m.sub,1/3,type=1))
		sup <- seq.int(quantile(block.l.m.sub,2/3,type=1),m)
		dif <- setdiff(block.l.m.sub ,c(inf,sup))
		block.l.m.sub <- c(sample(inf,1),if(!identical(dif,integer(0))) sample(dif,1),sample(sup,1))
								} 
	print(block.l.m.sub)
    cat(paste("Wait while the MBB resamples for the subsample are being computed;",type.sub.blocks,"computation enabled"))
	
	fit.mbb <- { 
		foreach(l=block.l.m.sub,.packages=c("gamlss.util","boot"),.errorhandling=c('pass'),.export=export)%:%foreach(j=m.subsamples)%op%{
		set.seed(seed) 
		switch(errorhandling
			,"try"= {
				MBB <- NULL
				attempt <- 1
				while( is.null(MBB) && attempt <= n.try ) {
				attempt <- attempt + 1
				try(MBB <- tsboot(subset(data,seq_len(n)%in%j), statistic=statistic, R = R, l = l, sim = 'fixed',...))
														   } 
					}
			,"pass"={MBB <- tsboot(subset(data,seq_len(n)%in%j), statistic=statistic, R = R, l = l, sim = 'fixed',...)})
	
		list(MBB=MBB,seed=seed,subsamble=j,
		index=seq_len(length(m.subsamples))[unlist(lapply(m.subsamples, function(subsamp) all(subsamp==j)))])                                 }
				}
	
	check01 <- lapply(seq_len(length(block.l.m.sub)), function(k)setdiff(seq_len(length(m.subsamples)),sapply(fit.mbb[[k]],function(j) j$index)) )
	if(!identical(unlist(check01),integer(0))) warning(sapply(seq_len(length(block.l.m.sub)), function(j) paste("length",block.l.m.sub[j],", iterations ",
	paste(as.character(if(identical(check01[[j]],integer(0)))0 else check01[[j]]),collapse=",")," have failed in step",iteration,";\n " ,sep="")))
	
	check <- lapply(seq_len(length(block.l.m.sub)), function(k) lapply(fit.mbb[[k]],function(j) if(is.null(j$MBB)) j$index))
	if(!identical(unlist(check),NULL)) warning(sapply(seq_len(length(block.l.m.sub)), function(j) paste("length",j,", iterations ",
	paste(sapply(check,unlist)[[j]],collapse=",")," have failed in step",iteration,";\n" ,sep="")))
	
	index1 <- block.l.m.sub[sapply(fit.mbb,length)!=0] ## for a given block length there might be a failure of the MBB
    
	sqe <- lapply(index1, function(k) {
	sqe <- sapply(seq_len(length(fit.mbb[[k]]))[setdiff(seq_len(length(fit.mbb[[k]])),unlist(check[[k]]))], function(j) {
	tn <- apply(fit.mbb[[k]][[j]]$MBB$t,2,mean)
	(t0-tn)^2})})

	(msqe <- lapply(sqe,function(j) {apply(j,1,mean)}))
	
	indice.best <- sapply(seq_len(length(msqe)), function(j) {mean(msqe[[j]][if(type.optm==0) TRUE else type.optm ])})
	l.m <- index1[which.min(indice.best)]
	
	opt.l$l.optm[iteration+1] <- round((n/length(block.l.m.sub))^(1/3)*l.m)
    
	cat(paste("Step:",iteration, "completed, more", nsteps-iteration,"to go"))
}
return(opt.l)
}
