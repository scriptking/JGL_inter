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

  ##############To identify C1 and C2 non-overlapping partition##############################
  
  if(K==2)  #use bi-conditional screening rule to identify block structure exactly
  {
    crit1 = list()
    for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 + lam2 }  
    S.sum = matrix(0,sum(connected),sum(connected))
    for(k in 1:K) {S.sum = S.sum + weights[k]*S[[k]]}
    S.sum = abs(S.sum)
    crit2 = S.sum > 2*lam1
  }
  
  if(K>2)  #use sufficient screening rule to identify larger-grained block structure
  {
    crit1 = list()
    for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 }  
    crit2 = matrix(0,sum(connected),sum(connected))
  }
  
  # are both criteria met?
  critboth = crit2
  for(k in 1:K) {critboth = critboth + crit1[[k]]}
  critboth = (critboth!=0)				
  diag(critboth) = 1

  
  #####################Create 2 nonoverlapping partitions C1 and C2############################
  
  ## now identify block structure using igraph:
  g1 <- graph.adjacency(critboth)	
  cout = clusters(g1)
  blocklist = list()
  # identify unconnected elements, and get blocks:
  unconnected = c()
  #	for(i in 2:(cout$no+1))
  #	{
  #		if(sum(cout$membership==(i-1))==1) { unconnected <- c(unconnected,which(cout$membership==(i-1))) }
  #		if(sum(cout$membership==(i-1))>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==(i-1)) }
  #	}
  
  # adapt cout$membership to start with index 1:
  if(min(cout$membership)==0){cout$membership=cout$membership+1}
  
  
  
  for(i in 1:(cout$no))
  {
    if(sum(cout$membership==i)==1) { unconnected <- c(unconnected,which(cout$membership==i)) }
    if(sum(cout$membership==i)>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==i) }
  }
  
  # final set of connected nodes
  connected[unconnected] = FALSE
  # connected indices of connected nodes:  0 for unconnected nodes, and 1:length(connected) for the rest.  
  # maps features 1:p to their order in the connected features
  connected.index = rep(0,p)
  connected.index[connected] = 1:sum(connected)
  # regular indices of connected nodes: map connected nodes onto 1:p indexing:
  
  # redefine unconnected as !connected (up until now it's been extra nodes caught as unconnected)
  unconnected=!connected
  
  ## define theta on all connected:   (so theta is really theta.connected).
  theta = list()
  for(k in 1:K) 
  {
    theta[[k]] = matrix(0,sum(connected),sum(connected))
    if(sum(connected)>0)
    {
      dimnames(theta[[k]])[[1]]=dimnames(theta[[k]])[[2]]=dimnames(Y[[k]])[[2]][connected]	
    }
  }
  
  ###################################################################################################
  
  
  
}








