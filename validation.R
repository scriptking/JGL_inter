validation <- function(inv_mat, valid, valid_type, ratio=0.1, inter_ratio=0.1,pred=0) {
  
  ### g1 means predict values of class 2 given values from values from class 1 i.e. conditioned on class 1 values
  ## mu_1g2 = sigma_12*sigma_22_inv*x_2
  ## mu_2g1 = sigma_21*sigma_11_inv*x_1
  
  p = dim(valid[[1]])[2]
  len_predict = floor(p*ratio)
  if(pred==0)
    to_predict = sample(1:p, len_predict, replace=F) ## Predict these features from rest
  else
    to_predict =pred
  len_predict = length(to_predict)
  
  if(valid_type==1) #intra
  {
    K = length(valid) # no of class
    
    err_percent = array(0,dim=c(K,len_predict))
    r2_value = array(0,dim=c(K,len_predict))
    
    for (k in 1:K) {
      val = (k - 1) * p
      sig = array(0, dim = c(p, p))
      
      i=1
      j=1
      for (i in 1:p) {
        for (j in 1:p) {
          sig[i, j] = inv_mat[val + i, val + j]
        }
      }
      
      temp_predict = array(0,dim=c(dim(valid[[k]])[1],len_predict))
      temp_ground = array(0,dim=c(dim(valid[[k]])[1],len_predict))
      temp2 = array(0,dim=c(1,len_predict))
      
      #to_predict = sample(1:p, len_predict, replace=F) ## Predict these features from rest
      ordered_sig <- sig[c(to_predict, setdiff(1:p, to_predict)), c(to_predict, setdiff(1:p, to_predict))] ## reorder the columns
      
      sigma_22_inv = solve(ordered_sig[ (len_predict + 1):p, (len_predict + 1):p ])
      sigma_12 = ordered_sig[1:len_predict, (len_predict + 1):p]
      
      temp_mul = sigma_12 %*% sigma_22_inv
      
      i=1
      j=1
      
      for (i in 1:dim(valid[[k]])[1]) {   # sample size in class
        
        pred_1g2 = temp_mul  %*% t(t(valid[[k]][i, setdiff(1:p, to_predict)]))
        ground_1g2 = t(valid[[k]][i, to_predict])
        
        temp2 = temp2 + (abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2)
        
        for(j in 1:len_predict){
          temp_predict[i,j] = pred_1g2[j,1]
          temp_ground[i,j] = ground_1g2[1,j]
        }
      }
      
      r2_val = r2error(temp_predict,temp_ground)
      temp2 = temp2 / i
      
      j=1
      for(j in 1:len_predict){
        err_percent[k,j] = temp2[1,j]
        r2_value[k,j] = r2_val[j]
      }
    }
    return(list(err_percent,r2_value))
  }
  
  if(valid_type==2) #inter
  {
    K = length(valid) # no of class
    total =  length(valid)*(length(valid)-1)/2 ##no of class combinations

    err_percent = array(0,dim=c(total,len_predict))
    r2_value = array(0,dim=c(total,len_predict))
    count = 0
    
    for (k1 in 1:(K-1)) {
      for (k2 in setdiff(k1:K,k1)) {
        val1 = (k1-1)*p
        val2 = (k2-1)*p
        count = count + 1
        
        sigma11 = array(0,dim=c(p,p))
        sigma22 = array(0,dim=c(p,p))
        sigma12 = array(0,dim=c(p,p))
        sigma = array(0,dim=c(2*p,2*p))
        
        i=1
        j=1
        for (i in 1:p) {
          for (j in 1:p) {
            sigma12[i,j] = inv_mat[val1+i,val2+j]
            sigma11[i,j] = inv_mat[val1+i,val1+j]
            sigma22[i,j] = inv_mat[val2+i,val2+j]
          }
        }
        i=1
        j=1
        for (i in 1:p) {
          for (j in 1:p) {
            sigma[i,j] = sigma11[i,j]
            sigma[p+i,p+j] = sigma22[i,j]
            sigma[i,p+j] = sigma12[i,j]
            sigma[p+i,j] = sigma12[j,i]
          }
        }
        
        temp_predict = array(0,dim=c(dim(valid[[k1]])[1],len_predict))
        temp_ground = array(0,dim=c(dim(valid[[k1]])[1],len_predict))
        temp2 = array(0,dim=c(1,len_predict))
        
        #to_predict = sample(1:p, len_predict, replace=F) ## Predict these features from rest
        ordered_sigma <- sigma[c(to_predict, setdiff(1:(2*p), to_predict)), c(to_predict, setdiff(1:(2*p), to_predict))] ## reorder the columns
        
        sigma_22_inv = solve(ordered_sigma[ (len_predict + 1):(2*p), (len_predict + 1):(2*p) ])
        sigma_12 = ordered_sigma[1:len_predict, (len_predict + 1):(2*p)]
        
        temp_mul = sigma_12 %*% sigma_22_inv
        
        i=1
        j=1
        
        for (i in 1:dim(valid[[k1]])[1]) {   # sample size in class
          
          pred_1g2 = temp_mul  %*% t(cbind(valid[[k1]][i, setdiff(1:p, to_predict)],valid[[k2]][i,1:p]))
          ground_1g2 = valid[[k1]][i, to_predict]
          
          temp2 = temp2 + (abs(t(pred_1g2) - ground_1g2)*100)/abs(ground_1g2)
          
          for(j in 1:len_predict){
            temp_predict[i,j] = pred_1g2[j,1]
            temp_ground[i,j] = ground_1g2[1,j]
          }
        }
        
        r2_val = r2error(temp_predict,temp_ground)
        temp2 = temp2 / i
        
        j=1
        for(j in 1:len_predict){
          err_percent[count,j] = temp2[1,j]
          r2_value[count,j] = r2_val[j]
        }
      }
    }
    
    return(list(err_percent,r2_value))
  }
}


## M2 = c(5,6,3,4,18,6,10,1,9,16,3,1,13,8,1,4,9,8,2,3,18,16,1,3,20)
## M3 = c(5,3)
## M1 = M1[c(M3,setdiff(1:5,M3)),c(M3,setdiff(1:5,M3))]