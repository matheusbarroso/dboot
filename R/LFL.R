LFL <-
function(data,statistic,R=100,iter=FALSE,nsteps=10,type.param=0,type.est="bias.variance",parallel=FALSE,...) {
if(is.logical(iter)==FALSE) stop("iter must be logical ") else if (is.numeric(nsteps)&(nsteps >0 )) nsteps <- ceiling(nsteps) else stop("nsteps must be a positive integer ") 

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
} else stop("unrecognised value of 'type.est',accepted values are: 'bias.variance' and 'distribution.quantile'")

n <- nrow(data)
c1 <- 1
R <- R

if (!iter) {
##step 1

l <- c1*n^(1/(r+4)) # or an iterative version
l <- ceiling(l)

estimate.l <- tsboot2(data, statistic=statistic, R = R, l = l, sim = 'fixed',...) ##paralelizar
estimate.2l <- tsboot2(data, statistic=statistic, R = R, l = 2*l, sim = 'fixed',...) ##parale

##step 2
endpt <- n - l + 1  ## endpoint
m <- c2*n^(1/3)*l^(2/3)  ## m= number of bootstrap blocks to be deleted
m <- ceiling (m)
if(m < endpt) "allgood" else "invalid value of m, m < n-3+1"
M <- endpt -m +1  
K <- 1:R
group <- M.collection(m=m,M=M,R=R,blocks=estimate.l$blocks) ## parallelizar

block.jab.point.value <- lapply(group, function(k) {

	if (is.null(k)) block.jab.point.value <- NULL else {
		x <- estimate.l$blocks$starts[k,]
		if(is.null(dim(x))) {
			blocks <- sapply(seq_len(dim(estimate.l$blocks$starts)[1]), function(j) {all(x%in%estimate.l$blocks$starts[j,])})              
			block.jab.point.value <-  estimate.l$t[blocks,]			
							} else {
			blocks <- apply(x,1,function (i) {sapply(seq_len(dim(estimate.l$blocks$starts)[1]), function(j) {all(i%in%estimate.l$blocks$starts[j,])})})
			boot.blocks <- apply(blocks,2,function(j) estimate.l$t[j,])
			block.jab.point.value <- apply(boot.blocks,1,mean)

	                               }
	                                                  }
})

block.jab.point.value <-matrix(unlist(block.jab.point.value),ncol=dim(estimate.l$t)[2],byrow=T)
pseudo <- (matrix(estimate.l$t0,ncol=dim(block.jab.point.value)[2],nrow=dim(block.jab.point.value)[1],byrow=T)*endpt-(endpt-m)*block.jab.point.value)/m 
var.jab <- apply(apply(pseudo,1,function(i)(i-estimate.l$t0)^2 ),1,sum)*m/((endpt-m)*M)

##step 3
C1.hat <- n*l^(-r)*var.jab
C2.hat <- 2*l*(estimate.l$t0-estimate.2l$t0)
l.optm <- (2*C2.hat^2/(r*C1.hat))^(1/(r+2))*n^(1/(r+2))
return(l.optm)} else {

l.init <- c1*n^(1/(r+4)) # or an iterative version
l.init <- ceiling(l.init)

length.path <- list(iteration=seq_len(nsteps),l=list())
length.path$l[[1]] <- l.init


for (iteration in seq_len(nsteps-1)) {

	##step 1

	l <- if(type.param==0) ceiling(mean(length.path$l[[iteration]])) else {
	if(is.numeric(type.param)&(type.param > 0)) length.path$l[[iteration]][type.param] else stop("'type.param' must be a positive integer,lesser than or equal to the number of parameters in the model.")
																		  }	

	estimate.l <- tsboot2(data, statistic=statistic, R = R, l = l, sim = 'fixed') ##paralelizar
	estimate.2l <- tsboot2(data, statistic=statistic, R = R, l = 2*l, sim = 'fixed') ##parale

	##step 2
	endpt <- n - l + 1  ## endpoint
	m <- c2*n^(1/3)*l^(2/3)  ## m= number of bootstrap blocks to be deleted
	m <- ceiling (m)
	if(m < endpt) "allgood" else "invalid value of m, m < n-3+1"
	M <- endpt -m +1  
	K <- 1:R
	group <- M.collection(m=m,M=M,R=R,blocks=estimate.l$blocks) ## parallelizar

	block.jab.point.value <- lapply(group, function(k) {

		if (is.null(k)) block.jab.point.value <- NULL else {
			x <- estimate.l$blocks$starts[k,]
			if(is.null(dim(x))) {
				blocks <- sapply(seq_len(dim(estimate.l$blocks$starts)[1]), function(j) {all(x%in%estimate.l$blocks$starts[j,])})              
				block.jab.point.value <-  estimate.l$t[blocks,]			
								} else {
				blocks <- apply(x,1,function (i) {sapply(seq_len(dim(estimate.l$blocks$starts)[1]), function(j) {all(i%in%estimate.l$blocks$starts[j,])})})
				boot.blocks <- apply(blocks,2,function(j) estimate.l$t[j,])
				block.jab.point.value <- apply(boot.blocks,1,mean)

										}
															}
														})

	block.jab.point.value <-matrix(unlist(block.jab.point.value),ncol=dim(estimate.l$t)[2],byrow=T)
	pseudo <- (matrix(estimate.l$t0,ncol=dim(block.jab.point.value)[2],nrow=dim(block.jab.point.value)[1],byrow=T)*endpt-(endpt-m)*block.jab.point.value)/m 
	var.jab <- apply(apply(pseudo,1,function(i)(i-estimate.l$t0)^2 ),1,sum)*m/((endpt-m)*M)

	##step 3
	C1.hat <- n*l^(-r)*var.jab
	C2.hat <- 2*l*(estimate.l$t0-estimate.2l$t0)
	l.optm <- (2*C2.hat^2/(r*C1.hat))^(1/(r+2))*n^(1/(r+2))
	length.path$l[[iteration+1]] <- l.optm^((r+2)/(r+4))
}
return(length.path)
}


}
