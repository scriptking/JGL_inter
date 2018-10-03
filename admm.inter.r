### ADMM for FGL:
admm.inter = function(S,lambda1,lambda2,rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5)
{
  K = length(S)
  p = dim(S[[1]])[2]
  #n=weights
  
  # initialize theta:
  theta = list()
  for(k in 1:K){theta[[k]] = 1+(0*(S[[k]]))}
  # initialize Z:
  Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
  # initialize W:
  W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }
  
  
  # iterations:
  iter=0
  diff_value = 10
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # reporting
    #	if(iter%%10==0)
    if(FALSE)
    {
      print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
      print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
    }
    
    
    # update theta:
    theta.prev = theta
    for(k in 1:K){
      edecomp = eigen(S[[k]] - rho*Z[[k]]/weights[k] + rho*W[[k]]/weights[k])
      D = edecomp$values
      V = edecomp$vectors
      D2 = weights[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/weights[k]) )
      theta[[k]] = Re(V %*% diag(D2) %*% t(V))
    }
    
    # update Z:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
    
    # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
    Z.prev = Z
    Z = flsa.inter.general(A,rho,lambda1,lambda2,penalize.diagonal=FALSE) # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
    
    # update the dual variable W:
    for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}
    
    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    diff_valuez = 0
    for(k in 1:K) {diff_valuez = diff_valuez + sum(abs(Z[[k]] - Z.prev[[k]])) / sum(abs(Z.prev[[k]]))}
    loss = 0; for(k in 1:K){loss = loss + sum(abs(theta[[k]]-Z[[k]]))}
    # increment rho by a constant factor:
    rho = rho*rho.increment
    print(sprintf("iter = %d,change_theta = %0.15f,change_z = %0.15f,loss = %0.15f",iter,diff_value,diff_valuez,loss))
  }
  
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  return(out)
}


