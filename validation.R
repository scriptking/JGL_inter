validation <- function(covar, valid) {
  inv_mat = solve(covar)
  intra_error = array(0, dim = length(valid))
  p = dim(valid[[1]])[2]
  intra_err = array(0,dim=c(length(valid),floor(p / 2)))
  for (k in 1:length(valid)) {
    val = (k - 1) * p
    sig = array(0, dim = c(p, p))
    for (i in 1:p) {
      for (j in 1:p) {
        sig[i, j] = inv_mat[val + i, val + j]
      }
    }
    
    sigma_22_1 = solve(sig[(floor(p / 2) + 1):p, (floor(p / 2) + 1):p])
    sigma_12 = sig[1:floor(p / 2), (floor(p / 2) + 1):p]
    val = 0
    temp = sigma_12 %*% sigma_22_1
    for (i in 1:dim(valid[[k]])[1]) {
      pred =  temp %*% t(valid[[k]][i, (floor(p / 2) + 1):p])
      ground = t(valid[[k]][i, 1:floor(p / 2)])
      val = val + sum(abs(ground - pred))
      intra_err[k,] = intra_err[k,] + abs(t(pred) - t(ground))
    }
    intra_error[k] = val / i
  }
  intra_err = intra_err / i
  
  K = length(valid)*(length(valid)-1)/2
  inter_error = array(0, dim = K)
  inter_err = array(0,dim=c(K,p))
  count = 0
  for (k1 in 1:(length(valid)-1)) {
    for (k2 in setdiff(k1:length(valid),k1)) {
        val1 = (k1-1)*p
        val2 = (k2-1)*p
        count = count + 1
        #sigma11 = array(0,dim=c(p,p))
        sigma22 = array(0,dim=c(p,p))
        sigma12 = array(0,dim=c(p,p))
        for (i in 1:p) {
          for (j in 1:p) {
            sigma12[i,j] = inv_mat[val1+i,val2+j]
            #sigma11[i,j] = inv_mat[((k1-1)*p)+i,((k1-1)*p)+j]
            sigma22[i,j] = inv_mat[((k2-1)*p)+i,((k2-1)*p)+j]
          }
        }

    sigma_22_1 = solve(sigma22)
    val = 0
    temp = sigma12 %*% sigma_22_1
    for (i in 1:dim(valid[[k1]])[1]) {
      pred =  temp %*% t(valid[[k2]][i,])
      ground = t(valid[[k1]][i, ])
      val = val + sum(abs(ground - pred))
      inter_err[count,] = inter_err[count,] + abs(t(pred) - t(ground))
    }
    inter_error[count] = val / i
    }
  }
  inter_err = inter_err / i
  return(list(intra_error,inter_error,intra_err,inter_err))
}