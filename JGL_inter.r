JGL_inter <- function(Y,lambda1=1,lambda2=1,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=500,tol=1e-5) {
  
  library(igraph)
  lambda1=1 
  lambda2=1 
  rho=1 
  weights="equal" 
  penalize.diagonal=FALSE 
  maxiter=500 
  tol=1e-5
  penalty="fused"
  
  
  load("~/JGL/example.data.rda")
  Y = example.data
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
  if(length(weights)==1){if(weights == "equal"){
    weights = rep(1,K)
  }}
  if(length(weights)==1){if(weights == "sample.size"){
    weights = rep(n/sum(n),K)
  }}
  
  connected = rep(TRUE,p)
  
  ## define S
  S = vector("list",length=K)
  for(k in 1:K)
  {
    ntemp = dim(Y[[k]])[1]
    #print(ntemp)
    S[[k]] = cov(Y[[k]][,connected])*(ntemp-1)/ntemp
  }
  
  
  # initialize lambdas:
  lam1 = penalty.as.matrix(lambda1,dim(Y[[1]])[2],penalize.diagonal=penalize.diagonal)
  lam2 = penalty.as.matrix(lambda2,dim(Y[[1]])[2],penalize.diagonal=TRUE)
  
  # run JGL on the block:
  Thetabl = admm.intra(S,lam1,lam2,penalty=penalty,rho=rho,weights=weights,penalize.diagonal=TRUE,maxiter=maxiter,tol=tol)
  # update Theta with Thetabl's results:
  for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}   
  
  
  
  
}








