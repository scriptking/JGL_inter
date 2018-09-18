library(igraph)
lambda1=1 
lambda2=1 
rho=1 
penalize.diagonal=FALSE 
maxiter=500 
tol=1e-5
load("~/JGL/example.data.rda")
Y = example.data

DATA = list()
DATA[[1]] = read.csv("~/JGL_inter/gene_adip.csv") ######1st column is identifier
DATA[[2]] = read.csv("~/JGL_inter/gene_musc.csv")
DATA[[3]] = read.csv("~/JGL_inter/gene_skin.csv")
DATA[[4]] = read.csv("~/JGL_inter/gene_thyr.csv")

K = 4
p_o = dim(DATA[[1]])[2]

vardata = array(dim=c(K,p_o-1))
for(k in 1:K)
{
  for(j in 2:p_o)
  {
    vardata[k,j-1] = var(DATA[[k]][,j])
  }
}
meandata = array(dim=c(1,p_o-1))
for(k in 1:K)
{
  for(j in 1:p_o-1)
  {
    meandata[j] = mean(vardata[,j])
  }
}

dim_name = dimnames(DATA[[1]])[[2]][2:p_o]
sorted_mean = sort(meandata,decreasing = TRUE, index.return=TRUE)
sorted_index = sorted_mean$ix
#DATA[[1]][2,sorted_index[1:10]+1]

DATA1 = list()
DATA1[[1]] = DATA[[1]][,sorted_index[1:100]+1]
DATA1[[2]] = DATA[[2]][,sorted_index[1:100]+1]
DATA1[[3]] = DATA[[3]][,sorted_index[1:100]+1]
DATA1[[4]] = DATA[[4]][,sorted_index[1:100]+1]
data1_dimname = dim_name[sorted_index[1:100]]

p = dim(DATA1[[1]])[2]

source('~/JGL_inter/JGL_inter.r')
source('~/JGL_inter/penalty.as.matrix.R')
source('~/JGL_inter/flsa.general.R')
source('~/JGL_inter/soft.R')
source('~/JGL_inter/admm.intra.r')

Z = JGL_inter(Y=DATA1)
print(Z$diff)