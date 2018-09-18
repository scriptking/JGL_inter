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
DATA[[1]] = read.csv("~/JGL_inter/gene_adip.csv")[,2:100]
DATA[[2]] = read.csv("~/JGL_inter/gene_musc.csv")[,2:100]
DATA[[3]] = read.csv("~/JGL_inter/gene_skin.csv")[,2:100]
DATA[[4]] = read.csv("~/JGL_inter/gene_thyr.csv")[,2:100]
K = 4
p = dim(DATA[[1]])[2]

vardata = array(dim=c(K,p))
for(k in 1:K)
{
  for(j in 1:p)
  {
    
    #vardata[k,j] = var(DATA[[k]][,j])
  }
}

source('~/JGL_inter/JGL_inter.r')
source('~/JGL_inter/penalty.as.matrix.R')
source('~/JGL_inter/flsa.general.R')
source('~/JGL_inter/soft.R')
source('~/JGL_inter/flsa.general.R')

Z = JGL_inter(DATA)
