JGL_inter <- function(Y,lambda1=1,lambda2=1,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=500,tol=1e-5) {
  
  library(igraph)
  lambda1=1 
  lambda2=1 
  rho=1 
  weights="sample.size" 
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
    weights = n/sum(n)
  }}
  
  connected = rep(TRUE,p)
  
  ## define S
  S = vector("list",length=K)
  for(k in 1:K)
  {
    ntemp = dim(Y[[k]])[1]
    print(ntemp)
    S[[k]] = cov(Y[[k]][,connected])*(ntemp-1)/ntemp
  }
  
  # if a penalty matrix is entered, only take its appropriate rows:
  lam1 = lambda1
  lam2 = lambda2

  

  for(i in 1:length(blocklist)){
      # the variables in the block
      bl <- blocklist[[i]] 
      Ybl = list()
      # get the data on only those variables
      for(k in 1:K) 
      {
        Ybl[[k]] = Y[[k]][,bl]
      }  
      # penalty matrices:
      if(length(lambda1)==1) { lam1.bl = lambda1 }
      if(length(lambda1)>1) { lam1.bl = lambda1[bl,bl] }
      if(length(lambda2)==1) { lam2.bl = lambda2 }
      if(length(lambda2)>1) { lam2.bl = lambda2[bl,bl] }
      # initialize lambdas:
      lam1.bl = penalty.as.matrix(lam1.bl,dim(Ybl[[1]])[2],penalize.diagonal=penalize.diagonal)
      if(penalty=="fused") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=TRUE)}
      if(penalty=="group") {lam2.bl = penalty.as.matrix(lam2.bl,dim(Ybl[[1]])[2],penalize.diagonal=penalize.diagonal)}
      
      # implement warm start if desired
      if(length(warm)==0) {warm.bl = NULL}
      if(length(warm)>0)
      {
        warm.bl = list()
        for(k in 1:K) { warm.bl[[k]] = warm[[k]][bl,bl] }
      }
      # run JGL on the block:
      Thetabl = admm.iters(Ybl,lam1.bl,lam2.bl,penalty=penalty,rho=rho,weights=weights,penalize.diagonal=TRUE,maxiter=maxiter,tol=tol,warm=warm.bl)
      # update Theta with Thetabl's results:
      for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}   
    }
  
  
  
}








