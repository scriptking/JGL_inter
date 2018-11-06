library(writexl)
library(corrplot)
source('~/JGL_inter/JGL_inter.r')
source('~/JGL_inter/penalty.as.matrix.R')
source('~/JGL_inter/flsa.general.R')
source('~/JGL_inter/soft.R')
source('~/JGL_inter/admm.intra.r')
source('~/JGL_inter/admm.inter.r')
source('~/JGL_inter/admm.iters.r')
source('~/JGL_inter/validation_intra.R')
source('~/JGL_inter/validation_inter.R')
source('~/JGL_inter/penalty.as.matrix.inter.R')
source('~/JGL_inter/flsa.inter.general.r')
source('~/JGL_inter/make.adj.matrix.R')
source('~/JGL_inter/net.degree.R')
source('~/JGL_inter/r2error.R')
source('~/JGL_inter/validation.R')

DATA = list()
DATA[[1]] = read.csv("~/JGL_inter/gene_adip.csv") ######1st column is identifier
DATA[[2]] = read.csv("~/JGL_inter/gene_musc.csv")
DATA[[3]] = read.csv("~/JGL_inter/gene_skin.csv")
DATA[[4]] = read.csv("~/JGL_inter/gene_thyr.csv")

tissues = array(c("adip","musc","skin","thyr"))

K = length(DATA)
p_o = dim(DATA[[1]])[2]

inter_issue = array(dim = (K*(K-1))/2)
count = 0
for (k1 in 1:(K-1)) {
  for (k2 in setdiff(k1:K,k1)) {
    count = count +1
    inter_issue[count] = paste(tissues[k1],tissues[k2],sep = "-")   
  }
}

for (k in 1:K) 
  DATA[[k]] = DATA[[k]][,2:p_o]

p = p_o - 1
vardata = array(dim=c(K,p))

for(k in 1:K)
  for(j in 1:p)
    vardata[k,j] = var(DATA[[k]][,j])

meandata = array(dim=c(1,p))
for(k in 1:K)
  for(j in 1:p)
    meandata[j] = mean(vardata[,j])

dim_name = dimnames(DATA[[1]])[[2]]
sorted_mean = sort(meandata,decreasing = TRUE, index.return=TRUE)
sorted_index = sorted_mean$ix

#for (k in 1:K) 
#  DATA[[k]] = DATA[[k]][,sorted_index[1:10000]]


DATA1 = list()
p = 40
for (k in 1:K)
  DATA1[[k]] = DATA[[k]][,sorted_index[1:p]]
data1_dimname = dim_name[sorted_index[1:p]]

l1.lineage = array(0,dim=c(1,(K*(K-1)/2)))
l1.lineage[1] = 2
l1.lineage[3] = 1
l1.lineage[5] = 1
l1.lineage[6] = 1
lambda1=3
lambda2=2
rho=1;penalize.diagonal=FALSE;maxiter=1000;tol=1e-5

DATA = list()

list_lambda = c(0.05,0.075,0.1,0.4,0.8,1.3,2,3)
#list_lambda = c(0.1,0.2)
count = 0

xl_col = c("iter","loss","l1","l2","total.edges.intra","total.edges.inter","avg.intra.unrelated.feature.percent.error.g2","avg.intra.unrelated.feature.r2.error.g2",
           "avg.intra.related.feature.percent.error.g2","avg.intra.related.feature.r2.error.g2")
observed_data_per = matrix(nrow = 1, ncol = length(xl_col), dimnames = list(1,xl_col))

for (x in 1:length(list_lambda)) {
  for (y in 1:length(list_lambda)) {
    count = count + 1
    lambda1=list_lambda[x]
    lambda2=list_lambda[y]
    Z = JGL_inter(Y=DATA1,l1_lineage=l1.lineage,lambda1,lambda2,rho,penalize.diagonal,maxiter,tol)
    validation.intra.nrelate = validation(Z$whole.sigma,Z$valid.data,1,0.2,pred=c(16,37,14,22))
    validation.intra.relate = validation(Z$whole.sigma,Z$valid.data,1,0.2,pred=c(10,40,39,7))
    degree.intra = net.degree(Z$theta.intra)
    degree.inter = net.degree(Z$theta.inter,inter=TRUE)
    #gra_intra = plot.jgl(Z$theta.intra)
    total.degree = list()
    total.degree[[1]] = rowSums(degree.intra)/2
    total.degree[[2]] = rowSums(degree.inter[[1]]+degree.inter[[2]])/2
    observed_data_per = rbind(observed_data_per,c(Z$iters,Z$diff,lambda1,lambda2,sum(total.degree[[1]]),sum(total.degree[[2]]),
                        mean(validation.intra.nrelate[[1]]),mean(validation.intra.nrelate[[2]]),
                        mean(validation.intra.relate[[1]]),mean(validation.intra.relate[[2]])))
    print(count)
    write_xlsx(data.frame(observed_data_per), "/home/manas/JGL_inter/valid_data_details.xlsx")
  }
}

#pred=c(10,40,39,34,24,7,26,30) ## corelated
#pred=c(2,12,16,37,14,22,4,8) ## uncorelated

observed_data_per = observed_data_per[2:dim(observed_data_per)[1],]
write_xlsx(data.frame(observed_data_per), "/home/manas/JGL_inter/valid_data_details.xlsx")

# corr=cor(DATA1[[1]])
# colnames(corr) <- 1:p
# out = corrplot(corr,order = "hclust", addrect = 8,hclust.method = "ward.D2")
# corrMatOrder(corr, order = "hclust", hclust.method = "ward.D2")
