validation_inter <- function(inv_mat, valid) {

  ### g1 means predict values of class 2 given values from values from class 1 i.e. conditioned on class 1 values

  ## mu_1g2 = sigma_12*sigma_22_inv*x_2
  ## mu_2g1 = sigma_21*sigma_11_inv*x_1
  
K = length(valid) # no of class
n =  length(valid)*(length(valid)-1)/2 ##no of class combinations
p = dim(valid[[1]])[2]

err_class = array(0, dim = n) #avg error in all features per inter class
err_g2 = array(0,dim=c(n,p)) ##avg error in all samples per inter class per feature [i,j] inter class i feature j
err_g1 = array(0,dim=c(n,p))
err_percent_g2 = array(0,dim=c(n,p))
err_percent_g1 = array(0,dim=c(n,p))
count = 0

for (k1 in 1:(K-1)) {
  for (k2 in setdiff(k1:K,k1)) {
    val1 = (k1-1)*p
    val2 = (k2-1)*p
    count = count + 1
    
    sigma11 = array(0,dim=c(p,p))
    sigma22 = array(0,dim=c(p,p))
    sigma12 = array(0,dim=c(p,p))
    
    for (i in 1:p) {
      for (j in 1:p) {
        sigma12[i,j] = inv_mat[val1+i,val2+j]
        sigma11[i,j] = inv_mat[val1+i,val1+j]
        sigma22[i,j] = inv_mat[val2+i,val2+j]
      }
    }
    
    sigma_22_inv = solve(sigma22)
    sigma_11_inv = solve(sigma11)
    val = 0
    
    temp_1g2 = sigma12 %*% sigma_22_inv
    temp_2g1 = t(sigma12) %*% sigma_11_inv
    
    temp1 = array(0,dim=c(1,p))
    temp2 = array(0,dim=c(1,p))
    temp3 = array(0,dim=c(1,p))
    temp4 = array(0,dim=c(1,p))
    
    for (i in 1:dim(valid[[k1]])[1]) { # sample size in class
      
      pred_1g2 =  temp_1g2 %*% t(valid[[k2]][i,])
      ground_1g2 = valid[[k1]][i,]
      
      #val = val + sum(abs(ground - pred))
      
      pred_2g1 =  temp_2g1 %*% t(valid[[k1]][i,])
      ground_2g1 = valid[[k2]][i,]
      
      temp1 = temp1 + abs(t(pred_1g2) - ground_1g2)
      temp2 = temp2 + abs(t(pred_2g1) - ground_2g1)
      temp3 = temp3 + ((abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2))
      temp4 = temp4 + ((abs(t(pred_2g1) - ground_2g1)*100)/abs(ground_2g1))
      
      #err_g2[count,] = err_g2[count,] + abs(t(pred_1g2) - ground_1g2)
      #err_g1[count,] = err_g1[count,] + abs(t(pred_2g1) - ground_2g1)
      
      #err_percent_g2[count,] = err_percent_g2[count,] + ((abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2))
      #err_percent_g1[count,] = err_percent_g1[count,] + ((abs(t(pred_2g1) - ground_2g1)*100)/abs(ground_2g1))
    }
    #inter_error[count] = val / i
    
    temp1 = temp1 / i
    temp2 = temp2 / i
    temp3 = temp3 / i
    temp4 = temp4 / i
    
    #err_g2[count,] = err_g2[count,] / i
    #err_g1[count,] = err_g1[count,] / i
    
    #err_percent_g2[count,] = err_percent_g2[count,] / i
    #err_percent_g1[count,] = err_percent_g1[count,] / i
    
    for(j in 1:p){
      err_g2[count,j] = temp1[1,j]
      err_g1[count,j] = temp2[1,j]
      err_percent_g2[count,j] = temp3[1,j]
      err_percent_g1[count,j] = temp4[1,j]
    }
  }
}

feature_err_g2 = colSums(err_g2)/dim(err_g2)[1] # avg error accross inter class for a feature
feature_err_g1 = colSums(err_g1)/dim(err_g1)[1]
feature_percent_err_g2 = colSums(err_percent_g2)/dim(err_percent_g2)[1] # avg percent error accross inter class for a feature
feature_percent_err_g1 = colSums(err_percent_g1)/dim(err_percent_g1)[1]

return(list(feature_percent_err_g2,feature_percent_err_g1))
}