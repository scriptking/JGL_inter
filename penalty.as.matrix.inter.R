penalty.as.matrix.inter <- function(lambda,p,K)
{
	# for matrix penalties:  check dim and symmetry:
	if(is.matrix(lambda))
	{
		if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
		if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
	}
	# for scalar penalties: convert to matrix form:
	if(length(lambda)==1) {
	  lambda=matrix(lambda,p,p)
	  lamda = vector("list",length=K)
	  for (i in 1:K)
	    lamda[[i]] = lambda
	  }
  
	return(lamda)
}

