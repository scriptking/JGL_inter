validation_intra <- function(inv_mat, valid) {
  
  ### g1 means predict values of class 2 given values from values from class 1 i.e. conditioned on class 1 values

  ## mu_1g2 = sigma_12*sigma_22_inv*x_2
  ## mu_2g1 = sigma_21*sigma_11_inv*x_1

  K = length(valid) # no of class
  p = dim(valid[[1]])[2]
  
  error = array(0, dim = K)
  err = array(0,dim=c(K,p))
  err_percent = array(0,dim=c(K,p))
  
  for (k in 1:K) {
    val = (k - 1) * p
    sig = array(0, dim = c(p, p))
    
    for (i in 1:p) {
      for (j in 1:p) {
        sig[i, j] = inv_mat[val + i, val + j]
      }
    }
    
    sigma_22_inv = solve(sig[ (floor(p / 2) + 1):p, (floor(p / 2) + 1):p ])
    sigma_11_inv = solve(sig[ 1:floor(p / 2), 1:floor(p / 2) ])
    sigma_12 = sig[1:floor(p / 2), (floor(p / 2) + 1):p]
    
    val = 0
    temp_g2 = sigma_12 %*% sigma_22_inv
    temp_g1 = t(sigma_12) %*% sigma_11_inv
    
    temp1 = array(0,dim=c(1,p))
    temp2 = array(0,dim=c(1,p))
    
    for (i in 1:dim(valid[[k]])[1]) {   # sample size in class

      pred_1g2 =  temp_g2 %*% t(valid[[k]][i, (floor(p / 2) + 1):p])
      ground_1g2 = valid[[k]][i, 1:floor(p / 2)]

      #val = val + sum(abs(ground - pred))

      pred_2g1 =  temp_g1 %*% t(valid[[k]][i, 1:floor(p / 2)])
      ground_2g1 = valid[[k]][i, (floor(p / 2) + 1):p]

      temp1 = temp1 + cbind(abs(t(pred_1g2) - ground_1g2),abs(t(pred_2g1) - ground_2g1))
      temp2 = temp2 + cbind((abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2),(abs(t(pred_2g1) - ground_2g1)*100)/abs(ground_2g1))
      #err[k,] = err[k,] + cbind(abs(t(pred_1g2) - ground_1g2),abs(t(pred_2g1) - ground_2g1))
      #err_percent[k,] = err_percent[k,] + cbind((abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2),(abs(t(pred_2g1) - ground_2g1)*100)/abs(ground_2g1))
    }
    #intra_error[k] = val / i
  #err[k,] = err[k,] / i
  #err_percent[k,] = err_percent[k,] / i
    temp1 = temp1 / i
    temp2 = temp2 / i
    
    for(j in 1:p){
      err[k,j] = temp1[1,j]
      err_percent[k,j] = temp2[1,j]
    }
  }

  feature_err = colSums(err)/dim(err)[1]
  feature_err_percent = colSums(err_percent)/dim(err_percent)[1]
  
  return(list(feature_err_percent))
}