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
