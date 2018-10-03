#load("~/JGL/example.data.rda")
#Y = example.data

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
source('~/JGL_inter/admm.inter.r')
source('~/JGL_inter/penalty.as.matrix.inter.R')
source('~/JGL_inter/flsa.inter.general.r')
source('~/JGL_inter/make.adj.matrix.R')
source('~/JGL_inter/net.degree.R')

l1.lineage = array(0,dim=c(1,(K*(K-1)/2)))
l1.lineage[2] = 1
l1.lineage[4] = 1
l1.lineage[6] = .5
lambda1=5
lambda2=1.5
rho=1.2
penalize.diagonal=FALSE 
maxiter=100 
tol=1e-5

Z = JGL_inter(Y=DATA1,l1_lineage=l1.lineage,lambda1,lambda2,rho,penalize.diagonal,maxiter,tol)
degree.intra = net.degree(Z$theta.intra)
degree.inter = net.degree(Z$theta.inter,inter=TRUE)
gra_intra = plot.jgl(Z$theta.intra)
count.disconnect = array(FALSE,dim=c(1,K+1))
for (k in 1:K) {
  count = 0
  for (i in 1:p) {
    if(!degree.intra[[k]][i])
      count = count + 1
  }
  count.disconnect[k] = count}
count = 0
for (k in 1:length(Z$theta.inter)) {
  for (i in 1:p) {
    if(!degree.inter[[k]][i])
      count = count + 1
  }}
count.disconnect[K+1] = count
count.disconnect

plot(gra_intra[[4]],layout=layout.fruchterman.reingold)#edge.width=E(graj[[4]])$weight*10)
adj = get.adjacency(gra_intra[[4]],attr='weight',sparse=FALSE)

plot(gra_intra[[4]],				#the graph to be plotted
     layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
     main='Organizational network example',	#specifies the title
     vertex.label.dist=0.5,			#puts the name labels slightly off the dots
     vertex.frame.color='blue', 		#the color of the border of the dots 
     vertex.label.color='black',		#the color of the name labels
     vertex.label.font=2,			#the font of the name labels
     vertex.label=NA, #V(graj[[4]])$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
     vertex.label.cex=1,			#specifies the size of the font of the labels. can also be made to vary
     vertex.size=5
)