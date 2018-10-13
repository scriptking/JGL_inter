net.degree <-
function(theta,inter=FALSE)
{
K = length(theta)
p = dim(theta[[1]])[1]
degree = array(0,dim=c(K,p))
for(k in 1:K)
{
  if(!inter){
	  degree[k,] = (rowSums(abs(theta[[k]])>1e-5)-1)
	  }
  else{
    degree[k,] = rowSums(abs(theta[[k]])>1e-5)
  }
}

return(degree)
}

