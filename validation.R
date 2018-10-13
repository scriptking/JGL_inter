validation <- function(covar, valid) {
  inv_mat = solve(covar)
  intra_error = array(0, dim = length(valid))
  p = dim(valid[[1]])[2]
  intra_err = array(0,dim=c(length(valid),p))
  for (k in 1:length(valid)) {
    val = (k - 1) * p
    sig = array(0, dim = c(p, p))
    for (i in 1:p) {
      for (j in 1:p) {
        sig[i, j] = inv_mat[val + i, val + j]
      }
    }
    
    sigma_22_1 = solve(sig[(floor(p / 2) + 1):p, (floor(p / 2) + 1):p])
    sigma_11_1 = solve(sig[1:floor(p / 2), 1:floor(p / 2)])
    sigma_12 = sig[1:floor(p / 2), (floor(p / 2) + 1):p]
    val = 0
    temp = sigma_12 %*% sigma_22_1
    temp1 = t(sigma_12) %*% sigma_11_1
    for (i in 1:dim(valid[[k]])[1]) {
      pred =  temp %*% t(valid[[k]][i, (floor(p / 2) + 1):p])
      ground = t(valid[[k]][i, 1:floor(p / 2)])
      val = val + sum(abs(ground - pred))
      pred1 =  temp1 %*% t(valid[[k]][i, 1:floor(p / 2)])
      ground1 = t(valid[[k]][i, (floor(p / 2) + 1):p])
      intra_err[k,] = intra_err[k,] + cbind(abs(t(pred) - t(ground)),abs(t(pred1) - t(ground1)))
    }
    intra_error[k] = val / i
  }
  intra_err = intra_err / i
  
  K = length(valid)*(length(valid)-1)/2
  inter_error = array(0, dim = K)
  inter_err = array(0,dim=c(K,p))
  inter_err1 = array(0,dim=c(K,p))
  count = 0
  for (k1 in 1:(length(valid)-1)) {
    for (k2 in setdiff(k1:length(valid),k1)) {
        val1 = (k1-1)*p
        val2 = (k2-1)*p
        count = count + 1
        sigma11 = array(0,dim=c(p,p))
        sigma22 = array(0,dim=c(p,p))
        sigma12 = array(0,dim=c(p,p))
        for (i in 1:p) {
          for (j in 1:p) {
            sigma12[i,j] = inv_mat[val1+i,val2+j]
            sigma11[i,j] = inv_mat[((k1-1)*p)+i,((k1-1)*p)+j]
            sigma22[i,j] = inv_mat[((k2-1)*p)+i,((k2-1)*p)+j]
          }
        }

    sigma_22_1 = solve(sigma22)
    sigma_11_1 = solve(sigma11)
    val = 0
    temp = sigma12 %*% sigma_22_1
    temp1 = t(sigma12) %*% sigma_11_1
    for (i in 1:dim(valid[[k1]])[1]) {
      pred =  temp %*% t(valid[[k2]][i,])
      ground = t(valid[[k1]][i, ])
      val = val + sum(abs(ground - pred))
      pred1 =  temp1 %*% t(valid[[k1]][i,])
      ground1 = t(valid[[k2]][i, ])
      inter_err[count,] = inter_err[count,] + abs(t(pred) - t(ground))
      inter_err1[count,] = inter_err1[count,] + abs(t(pred1) - t(ground1))
      
    }
    inter_error[count] = val / i
    }
  }
  inter_err = inter_err / i
  inter_err1 = inter_err1 / i
  intra_gene_err = colSums(intra_err)/dim(intra_err)[1]
  inter_gene_err = colSums(inter_err)/dim(inter_err)[1]
  inter_gene_err1 = colSums(inter_err1)/dim(inter_err1)[1]
  return(list(intra_error,inter_error,intra_gene_err,inter_gene_err,inter_gene_err1,intra_err,inter_err,inter_err1))
}