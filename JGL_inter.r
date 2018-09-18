JGL_inter <- function(Y,lambda1=1,lambda2=1,rho=1,penalize.diagonal=TRUE,maxiter=500,tol=1e-5) {
  
  p = dim(Y[[1]])[2] # No of genes  i.e. dimention size
  K = length(Y) # No of tissues i.e. total classes
  sample_size = dim(Y[[1]])[1]
  n = sample_size #total no of samples
  #n = rep(0,K)
  print(sprintf("dimention=%d, classes=%d, sample=%d", p,K,sample_size))
  #for(k in 1:K) {n[k] = dim(Y[[k]])[1]}---unequal sample size store in n[1] to n[K] variable
  
  # assign feature names if none exist:
  if(length(dimnames(Y[[1]])[[2]])==0)
  {
    for(k in 1:K)
    {
      dimnames(Y[[k]])[[2]]=paste("gene ",1:p,sep="")
    }
  }
  
  # mean-normalize Y:
  for(k in 1:K)
  {
    for(j in 1:p)
    {
      Y[[k]][,j] = Y[[k]][,j]-mean(Y[[k]][,j])
    }
  }
  
  # set weights to each class:
  weights = rep(n,K)
  
  ## define S_intra
  S_intra = vector("list",length=K)
  for(k in 1:K)
  {
    S_intra[[k]] = cov(Y[[k]])*(n-1)/n
  }
  
  
  ## define S_inter
  S_inter = vector("list",length=K)
  k = 0 
  S_inter_seq = array(dim=c((K*(K-1))/2,2))
  for(k1 in 1:K)
  {
    for (k2 in setdiff(k1:K,k1)) 
    {
      #print(sprintf("k1 = %d, k2 = %d",k1,k2))
      k = k + 1
      S_inter[[k]] = cov(Y[[k1]],Y[[k2]])*(n-1)/n
      S_inter_seq[k,1] = k1
      S_inter_seq[k,2] = k2
    }
  }
  
  
  # initialize lambdas:
  lam1 = penalty.as.matrix(lambda1,dim(Y[[1]])[2],penalize.diagonal=penalize.diagonal)
  lam2 = penalty.as.matrix(lambda2,dim(Y[[1]])[2],penalize.diagonal=TRUE)
  
  # run JGL for intra class:
  print("before admm call for intra")
  Theta = admm.intra(S_intra,lam1,lam2,rho=rho,penalize.diagonal=TRUE,weights=weights,maxiter=maxiter,tol=tol)
  # update Theta with Theta's results:
  theta_intra = list()
  for(k in 1:K) {theta_intra[[k]] = Theta$Z[[k]]}   
  
  
  return(theta_intra)
}








