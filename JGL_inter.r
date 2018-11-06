JGL_inter <- function(Y,l1_lineage,lambda1=1,lambda2=1,rho=1,penalize.diagonal=FALSE,maxiter=500,tol=1e-5,pred=0) {
  
  p = dim(Y[[1]])[2]              # No of genes  i.e. dimention size
  K = length(Y)                   # No of tissues i.e. total classes
  sample_size = dim(Y[[1]])[1]
  n = sample_size                 #total no of samples

  
  # assign feature names if none exist:
  if(length(dimnames(Y[[1]])[[2]])==0)
    for(k in 1:K)
      dimnames(Y[[k]])[[2]]=paste("dim_",1:p,sep="")
    
  # mean-normalize Y:
  for(k in 1:K)
    Y[[k]] = scale(Y[[k]], center = TRUE, scale = TRUE)
    #for(j in 1:p)
      #Y[[k]][,j] = Y[[k]][,j]-mean(Y[[k]][,j])
  
  
  
  Y_valid = vector("list",length=K)
  n_train = floor(n * .9)
  for(k in 1:K)
  {
      Y_valid[[k]] = Y[[k]][1+n_train:n-1,]
      Y[[k]] = Y[[k]][1:n_train,]
  }
  # set weights to each class:
  n = n_train
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
  for(k1 in 1:(K-1))
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
  Theta = admm.iters(S_intra,S_inter,lam1,lam2,rho=rho,penalize.diagonal=TRUE,weights=weights,maxiter=maxiter,tol=tol)
  
  # round very small theta entries down to zero:
  rounddown = abs(Theta$Z)<tol;
  Theta$Z=Theta$Z*(1-rounddown)
  
  theta_intra = vector("list",length=K)
  for (k in 1:K) {
    val = (k-1)*p
    theta_intra[[k]] = array(0,dim=c(p,p))
    for (i in 1:p) {
      for (j in 1:p) {
        theta_intra[[k]][i,j] = Theta$Z[val+i,val+j]
      }
    }
  }
  theta_inter = vector("list",length=K*(K-1)/2)
  count = 0
  for (k1 in 1:(K-1)) {
    for (k2 in setdiff(k1:K,k1)) {
      val1 = (k1-1)*p
      val2 = (k2-1)*p
      count = count + 1
      theta_inter[[count]] = array(0,dim=c(p,p))
      for (i in 1:p) {
        for (j in 1:p) {
          theta_inter[[count]][i,j] = Theta$Z[val1+i,val2+j]
        }
      }
    }
  }
  inverse_covar = solve(Theta$Z)
  #validation.total_inter = validation(inverse_covar,Y_valid,2,0.2)
  out = list(whole.sigma=inverse_covar,valid.data=Y_valid,theta.intra=theta_intra,diff=Theta$diff,iters=Theta$iters,theta.inter=theta_inter)
  
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
