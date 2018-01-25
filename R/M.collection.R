#'Internal \code{\link{LFL}} function for computing I(i)
#' in the JAB variance estimator.
#'
#'\code{M.collection} computes the subcollection of all
#'resampled block-sets of size \code{R} where none of the
#'resampled blocks equals the deleted blocks.
#'
#'\code{M.collection} is an internal function called by
#'\code{\link{LFL}}. It has the resampled blocks as inputs
#'and the reusable blocks as output.
#'
#'@author Matheus de Vasconcellos Barroso
#'
#
#'

#'@param m The number of bootstrap blocks
#' to be deleted. The default value is:
#' \eqn{m = c_{2}*n^(1/3)*l^(2/3)}
#'
#'@param M The number of usable observations after
#'deleting \code{m}, \eqn{M = n-l-m-2}
#'
#'@param R A positive integer giving the number of
#' bootstrap replicates required.
#'
#'@param blocks A list with two elements, the first named
#'starts and the second lenghts. The former is a pxq matrix,
#'\eqn{p = R} and \eqn{q = ceiling(n/l)} giving the p block
#'resample initial positions. The latter gives the length of
#'each block. \code{blocks} is usually equal to
#' \code{estimate.l$blocks}, the block output of the
#'  \code{\link{tsboot2}} function.
#'
#'@return A list with the index of reusable blocks of length
#' \code{M}. This is the I(i) step in the JAB algorithm.
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


M.collection <-
function(m,M,R,blocks) lapply(1:M, function(i) {
	Bi <- i:(i+m-1) ## generating the sequence B(i),...,B(i+m-1)
	Ii <- sapply(seq_len(R), function(k) {
		ends <-  cbind(blocks$starts[k, ], blocks$lengths) # start and length of each block in the resample
		check <- Bi%in%ends[,1] # check whether there is a matching element (should be block) from Bi in ends
		check <- TRUE%in%check #at least one Bi in ends
		if(!check) k  # if there is no Bi in ends this behaves as a r.s.
                                          })

	Ii <- unlist(Ii)

})
