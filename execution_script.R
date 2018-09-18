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
p = 10
for (k in 1:K) {
  DATA1[[k]] = DATA[[k]][,sorted_index[1:p]+1]
}
data1_dimname = dim_name[sorted_index[1:p]]

if(FALSE) {
DATA1[[1]] = DATA[[1]][,sorted_index[1:100]+1]
DATA1[[2]] = DATA[[2]][,sorted_index[1:100]+1]
DATA1[[3]] = DATA[[3]][,sorted_index[1:100]+1]
DATA1[[4]] = DATA[[4]][,sorted_index[1:100]+1]
data1_dimname = dim_name[sorted_index[1:100]]
p = dim(DATA1[[1]])[2]
}

source('~/JGL_inter/JGL_inter.r')
source('~/JGL_inter/penalty.as.matrix.R')
source('~/JGL_inter/flsa.general.R')
source('~/JGL_inter/soft.R')
source('~/JGL_inter/admm.intra.r')
source('~/JGL_inter/make.adj.matrix.R')

Z = JGL_inter(Y=DATA1)
theta = Z$theta
print(Z$diff)
plot(graj[[k]])

K=length(theta)
adj = make.adj.matrix(theta,separate=TRUE)
for (k in 1:K) {
  diag(adj[[k]])=0
}

graj = list()
for (k in 1:K) {
  graj[[k]]= graph.adjacency(adj[[k]],mode="upper",weighted=TRUE)
}
V(gadj)$label <- data1_dimname[V(gadj)]
#gadj = graph.adjacency(adj,mode="upper",weighted=TRUE)
#weight the edges according to the classes they belong to
E(gadj)$color = 2^(K)-get.edge.attribute(gadj,"weight")
#plot the net using igraph
plot(gadj, vertex.frame.color="white",layout=layout.fruchterman.reingold, 
     vertex.label=NA, vertex.label.cex=3, vertex.size=1)
#edge_weight = theta[[1]][i,j]
count = 1
edge_weight=array(dim = c(1,length(E(gadj))))
for (i in 1:(p-1)) {
  for (j in (i+1):p) {
    if(adj[[1]][i,j]){
          #print(sprintf("loop %d-- %d,,,,count = %d",i,j,count))
          #print(E(gadj)[count])
          edge_weight[count] = theta[[1]][i,j]
          count = count +1
      }
      
  }
}
#V(gadj)$label <- data1_dimname[V(gadj)]
V(gadj)$name <- data1_dimname[V(gadj)]
E(gadj)$weight = edge_weight
plot(gadj)