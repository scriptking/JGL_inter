#load("~/JGL/example.data.rda")
#Y = example.data

DATA = list()
DATA[[1]] = read.csv("~/JGL_inter/gene_adip.csv") ######1st column is identifier
DATA[[2]] = read.csv("~/JGL_inter/gene_musc.csv")
DATA[[3]] = read.csv("~/JGL_inter/gene_skin.csv")
DATA[[4]] = read.csv("~/JGL_inter/gene_thyr.csv")

tissues = array(c("adip","musc","skin","thyr"))

K = 4
p_o = dim(DATA[[1]])[2]

inter_issue = array(dim = (K*(K-1))/2)
count = 0
for (k1 in 1:(K-1)) {
  for (k2 in setdiff(k1:K,k1)) {
    count = count +1
    inter_issue[count] = paste(tissues[k1],tissues[k2],sep = "-")   
  }
}

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
source('~/JGL_inter/admm.iters.r')
source('~/JGL_inter/validation.R')
source('~/JGL_inter/penalty.as.matrix.inter.R')
source('~/JGL_inter/flsa.inter.general.r')
source('~/JGL_inter/make.adj.matrix.R')
source('~/JGL_inter/net.degree.R')

l1.lineage = array(0,dim=c(1,(K*(K-1)/2)))
l1.lineage[1] = 2
l1.lineage[3] = 1
l1.lineage[5] = 1
l1.lineage[6] = 1
lambda1=3
lambda2=2
rho=1;penalize.diagonal=FALSE;maxiter=1000;tol=1e-5

DATA = list()
#list_lambda = c(0.1,0.4,0.8,1.1,1.5,2,3,5)
list_lambda = c(0.1,0.2)
out_data = list()
count = 0
col_name = c("iter","loss","l1","l2",paste("deg.",tissues,sep=""),paste("deg.",inter_issue,sep=""),paste("intra.gene.err.g2.",data1_dimname[1:floor(p/2)],sep=""),paste("intra.gene.err.g1.",data1_dimname[(1+floor(p/2)):p],sep=""),paste("inter.gene.err.g2.",data1_dimname,sep=""),paste("inter.gene.err.g1.",data1_dimname,sep=""))
observed_data = matrix(nrow = 1, ncol = length(col_name), dimnames = list(1,col_name))
for (x in 1:length(list_lambda)) {
  for (y in 1:length(list_lambda)) {
    count = count + 1
    lambda1=list_lambda[x]
    lambda2=list_lambda[y]
    Z = JGL_inter(Y=DATA1,l1_lineage=l1.lineage,lambda1,lambda2,rho,penalize.diagonal,maxiter,tol)
    #print(sprintf("iter=%d,loss = %f,l1=%f,l2=%f",Z$iters,Z$diff,lambda1,lambda2))
    #print(Z$validation.error[[1]])
    #print(Z$validation.error[[2]])
    #print(Z$validation.error[[3]])
    #print(Z$validation.error[[4]])
    degree.intra = net.degree(Z$theta.intra)
    degree.inter = net.degree(Z$theta.inter,inter=TRUE)
    gra_intra = plot.jgl(Z$theta.intra)
    #count.disconnect = array(FALSE,dim=c(1,K+1))
    total.degree = list()
    total.degree[[1]] = rowSums(degree.intra)
    total.degree[[2]] = rowSums(degree.inter[[1]])
    #print(total.degree[[1]])
    #print(total.degree[[2]])
    out_str = paste(sprintf("iter=%d,loss=%f,l1=%f,l2=%f,",Z$iters,Z$diff,lambda1,lambda2),"intra.deg=",toString(total.degree[[1]]),"inter.deg=",toString(total.degree[[2]]),"intra.err=",toString(sprintf("%f",Z$validation.error[[1]])),"inter.err=",toString(sprintf("%f",Z$validation.error[[2]])))
    print(out_str)
    observed_data = rbind(observed_data,c(Z$iters,Z$diff,lambda1,lambda2,total.degree[[1]],total.degree[[2]],Z$validation.error[[3]],Z$validation.error[[4]],Z$validation.error[[5]]))
    #c("iter","loss","l1","l2",paste("deg.",tissues,sep=""),paste("deg.",inter_issue,sep=""),paste("val.err.",tissues,sep=""),paste("val.err.",inter_issue,sep=""),paste("intra.gene.err.",data1_dimname[1:length(Z$validation.error[[3]])],sep=""),paste("inter.gene.err.",data1_dimname[1:length(Z$validation.error[[4]])],sep=""))
  }
}
observed_data[2:dim(observed_data)[1],]

plot(gra_intra[[4]],layout=layout.fruchterman.reingold)#edge.width=E(graj[[4]])$weight*10)
adj = get.adjacency(gra_intra[[4]],attr='weight',sparse=FALSE)

plot(gra_intra[[4]],				#the graph to be plotted
     layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
     main='layer 4',	#specifies the title
     vertex.label.dist=0.5,			#puts the name labels slightly off the dots
     vertex.frame.color='blue', 		#the color of the border of the dots 
     vertex.label.color='black',		#the color of the name labels
     vertex.label.font=2,			#the font of the name labels
     vertex.label=NA, #V(graj[[4]])$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
     vertex.label.cex=1,			#specifies the size of the font of the labels. can also be made to vary
     vertex.size=5
)