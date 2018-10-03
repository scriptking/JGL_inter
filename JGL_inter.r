JGL_inter <- function(Y,l1_lineage,lambda1=1,lambda2=1,rho=1,penalize.diagonal=FALSE,maxiter=500,tol=1e-5) {
  
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
  S_inter = vector("list",length=(K*(K-1))/2)
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
  inter_lam1 = penalty.as.matrix.inter(lambda1,dim(Y[[1]])[2],length(S_inter))
  inter_lam2 = penalty.as.matrix.inter(lambda2,dim(Y[[1]])[2],length(S_inter))
  
  l1.lineage = vector("list",length=length(l1_lineage))
  for (i in 1:length(l1_lineage)) {
    l1.lineage[[i]] = penalty.as.matrix(l1_lineage[i],dim(Y[[1]])[2],penalize.diagonal=TRUE)
    inter_lam1[[i]] = inter_lam1[[i]] + l1.lineage[[i]]
  }
  
  # run JGL for intra class:
  print("before admm call for intra")
  Theta = admm.intra(S_intra,lam1,lam2,rho=rho,penalize.diagonal=TRUE,weights=weights,maxiter=maxiter,tol=tol)
  Theta_inter = admm.inter(S_inter,inter_lam1,inter_lam2,rho=rho,weights=rep(n,length(S_inter)),maxiter=maxiter,tol=tol)
  
  # update Theta with Theta's results:
  theta_intra = list()
  theta_inter = list()
  for(k in 1:K) {theta_intra[[k]] = Theta$Z[[k]]}   
  for(k in 1:length(S_inter)) {theta_inter[[k]] = Theta_inter$Z[[k]]} 
  
  # round very small theta entries down to zero:
  for(k in 1:K)
  {
    rounddown = abs(theta_intra[[k]])<tol; diag(rounddown)=FALSE
    theta_intra[[k]]=theta_intra[[k]]*(1-rounddown)
  }
  
  for(k in 1:length(S_inter))
  {
    rounddown = abs(theta_inter[[k]])<tol;
    theta_inter[[k]]=theta_inter[[k]]*(1-rounddown)
  }
  
  out = list(theta.intra=theta_intra,diff.intra=Theta$diff,iters.intra=Theta$iter,theta.inter=theta_inter,diff.inter=Theta_inter$diff,iters.inter=Theta_inter$iter)
  
  
  return(out)
}

plot.jgl <- function(x,...)
  {
    theta=x
    library(igraph)
    K=length(theta)
    adj = make.adj.matrix(theta,separate=TRUE)
    for (k in 1:K) { diag(adj[[k]])=0 }
    graj = list()
    for (k in 1:K) { graj[[k]]= graph.adjacency(adj[[k]],mode="upper",weighted=TRUE) }

    edge_weight = list()
    for (k in 1:K) { edge_weight[[k]]=array(dim = c(1,length(E(graj[[k]])))) }
    for (k in 1:K) {
      count = 1
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          if(adj[[k]][i,j]){
            edge_weight[[k]][count] = theta[[k]][i,j]
            count = count +1
          }}}}
    
    for (k in 1:K) {
    #V(graj[[k]])$name <- data1_dimname[V(graj[[k]])]
    E(graj[[k]])$weight = edge_weight[[k]]
    }
  
    return(graj)
}
