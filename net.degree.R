net.degree <-
function(theta,inter=FALSE)
{
K = length(theta)
degree = list()
for(k in 1:K)
{
  if(!inter){
	  degree[[k]] = (rowSums(abs(theta[[k]])>1e-5)-1)
	  }
  else{
    degree[[k]] = rowSums(abs(theta[[k]])>1e-5)
  }
}

return(degree)
}

