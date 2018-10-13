### ADMM for FGL:
admm.iters = function(S_intra,S_inter,lambda1,lambda2,rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5)
{
  K = length(S_intra)
  p = dim(S_intra[[1]])[2]
  n=weights
  
  S_big = array(0,dim=c(K*p,K*p))
  for (k in 1:K) {
    val = (k-1)*p
    for (i in 1:p) {
      for (j in 1:p) {
        S_big[val+i,val+j] = S_intra[[k]][i,j]
      }
    }
  }
  count = 0
  for (k1 in 1:(K-1)) {
    for (k2 in setdiff(k1:K,k1)) {
      val1 = (k1-1)*p
      val2 = (k2-1)*p
      count = count + 1
      for (i in 1:p) {
        for (j in 1:p) {
          S_big[val1+i,val2+j] = S_inter[[count]][i,j]
          S_big[val2+j,val1+i] = S_inter[[count]][i,j]
        }
      }
    }
  }
  # initialize theta:
  theta = array(1,dim=c(K*p,K*p))
  # initialize Z:
  Z = array(0,dim=c(K*p,K*p))
  # initialize W:
  W = array(0,dim=c(K*p,K*p))
  
  
  # iterations:
  iter=0
  diff_value = 10
  diff_valuez = 0
  
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # reporting
    #	if(iter%%10==0)
    if(FALSE)
    {
      #print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
      #print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
    }
    
    
    # update theta:
    theta.prev = theta
    edecomp = eigen(S_big - rho*Z/n + rho*W/n)
    D = edecomp$values
    V = edecomp$vectors
    D2 = n/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n) )
    theta = V %*% diag(D2) %*% t(V)
    
    
    # update Z:
    # define A matrices:
    A = theta + W
    
    # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
    Z.prev = Z
    
    A_intra = vector("list",length=K)
    for (k in 1:K) {
      val = (k-1)*p
      A_intra[[k]] = array(0,dim=c(p,p))
      for (i in 1:p) {
        for (j in 1:p) {
          A_intra[[k]][i,j] = A[val+i,val+j]
        }
      }
    }
    A_inter = vector("list",length=K*(K-1)/2)
    count = 0
    for (k1 in 1:(K-1)) {
      for (k2 in setdiff(k1:K,k1)) {
        val1 = (k1-1)*p
        val2 = (k2-1)*p
        count = count + 1
        A_inter[[count]] = array(0,dim=c(p,p))
        for (i in 1:p) {
          for (j in 1:p) {
            A_inter[[count]][i,j] = A[val1+i,val2+j]
          }
        }
      }
    }
    
    #if(K==2){Z = flsa2(A,rho,lambda1,lambda2,penalize.diagonal=TRUE)}
    Z_intra = flsa.general(A_intra,rho,lambda1,lambda2,penalize.diagonal=TRUE)  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
    Z_inter = flsa.inter.general(A_inter,rho,lambda1,lambda2,penalize.diagonal=FALSE) # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
    # update the dual variable W:
    
    Z = array(0,dim=c(K*p,K*p))
    for (k in 1:K) {
      val = (k-1)*p
      for (i in 1:p) {
        for (j in 1:p) {
          Z[val+i,val+j] = Z_intra[[k]][i,j]
        }
      }
    }
    count = 0
    for (k1 in 1:(K-1)) {
      for (k2 in setdiff(k1:K,k1)) {
        val1 = (k1-1)*p
        val2 = (k2-1)*p
        count = count + 1
        for (i in 1:p) {
          for (j in 1:p) {
            Z[val1+i,val2+j] = Z_inter[[count]][i,j]
            Z[val2+j,val1+i] = Z_inter[[count]][i,j]
          }
        }
      }
    }
    W = W + (theta-Z)
    
    # bookkeeping:
    iter = iter+1
    diff_value = sum(abs(theta - theta.prev)) / sum(abs(theta.prev))
    diff_valuez = sum(abs(Z - Z.prev)) / sum(abs(Z.prev))
    loss = sum(abs(theta - Z))
    # increment rho by a constant factor:
    rho = rho*rho.increment
    #print(sprintf("iter = %d,change_theta = %0.15f,change_z = %0.15f,loss = %0.15f",iter,diff_value,diff_valuez,loss))
  }
  
  diff = sum(abs(theta - Z))
  out = list(theta=theta,Z=Z,theta_diff=diff_value,Z_diff=diff_valuez,diff=loss,iters=iter)
  
  return(out)
}


