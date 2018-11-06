r2error <- function(prediction,ground){

rss <- colSums((prediction - ground) ^ 2)  ## residual sum of squares
tss <- colSums((ground - colMeans(ground)) ^ 2)  ## total sum of squares
rsq <- 1 - (rss/tss)

return(rsq)
}