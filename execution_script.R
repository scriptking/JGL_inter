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
p = 10
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

#list_lambda = c(0.05,0.075,0.1,0.4,0.8,1.3,2,3)
list_lambda = c(0.1,0.2)
out_data = list()
count = 0
#col_name = c("iter","loss","l1","l2",paste("deg.",tissues,sep=""),paste("deg.",inter_issue,sep=""),paste("intra.gene.err.g2.",data1_dimname[1:floor(p/2)],sep=""),paste("intra.gene.err.g1.",data1_dimname[(1+floor(p/2)):p],sep=""),paste("inter.gene.err.g2.",data1_dimname,sep=""),paste("inter.gene.err.g1.",data1_dimname,sep=""))
col_name = c("iter","loss","l1","l2",paste("deg.",tissues,sep=""),paste("tra.g2.",data1_dimname[1:floor(p/2)],sep=""),paste("tra.g1.",data1_dimname[(1+floor(p/2)):p],sep=""),paste("ter.g2.",data1_dimname,sep=""),paste("ter.g1.",data1_dimname,sep=""))
percent_err_col_name = c("iter","loss","l1","l2",paste("tra.g2.per.",data1_dimname[1:floor(p/2)],sep=""),paste("tra.g1.per.",data1_dimname[(1+floor(p/2)):p],sep=""),paste("ter.avg.per.",data1_dimname,sep=""),paste("ter.g2.per.",data1_dimname,sep=""),paste("ter.g1.per.",data1_dimname,sep=""))

observed_data = matrix(nrow = 1, ncol = length(col_name), dimnames = list(1,col_name))
observed_data_per = matrix(nrow = 1, ncol = length(percent_err_col_name), dimnames = list(1,percent_err_col_name))

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
    #out_str = paste(sprintf("iter=%d,loss=%f,l1=%f,l2=%f,",Z$iters,Z$diff,lambda1,lambda2),"intra.deg=",toString(total.degree[[1]]),"inter.deg=",toString(total.degree[[2]]),"intra.err=",toString(sprintf("%f",Z$validation.error[[1]])),"inter.err=",toString(sprintf("%f",Z$validation.error[[2]])))
    #print(out_str)
    #observed_data = rbind(observed_data,c(Z$iters,Z$diff,lambda1,lambda2,total.degree[[1]],total.degree[[2]],Z$validation.error[[3]],Z$validation.error[[4]],Z$validation.error[[5]]))
    #c("iter","loss","l1","l2",paste("deg.",tissues,sep=""),paste("deg.",inter_issue,sep=""),paste("val.err.",tissues,sep=""),paste("val.err.",inter_issue,sep=""),paste("intra.gene.err.",data1_dimname[1:length(Z$validation.error[[3]])],sep=""),paste("inter.gene.err.",data1_dimname[1:length(Z$validation.error[[4]])],sep=""))
    #observed_data = rbind(observed_data,c(Z$iters,Z$diff,lambda1,lambda2,total.degree[[1]],Z$validation.error[[3]],Z$validation.error[[4]],Z$validation.error[[5]]))
    avg_err_per_inter = (Z$validation.inter[[1]] + Z$validation.inter[[2]])/2
    observed_data_per = rbind(observed_data_per,c(Z$iters,Z$diff,lambda1,lambda2,Z$validation.intra[[1]],avg_err_per_inter,Z$validation.inter[[1]],Z$validation.inter[[2]]))
    print(count)
    #write_xlsx(data.frame(observed_data), "/home/manas/JGL_inter/observed_data_error.xlsx")
    #write_xlsx(data.frame(observed_data_per), "/home/manas/JGL_inter/observed_data_percent_1.xlsx")
  }
}
#observed_data = observed_data[2:dim(observed_data)[1],]
observed_data_per = observed_data_per[2:dim(observed_data_per)[1],]
#print(observed_data)

#write_xlsx(data.frame(observed_data), "/home/manas/JGL_inter/observed_data_error.xlsx")
write_xlsx(data.frame(observed_data_per), "/home/manas/JGL_inter/observed_data_percent_1.xlsx")

#plot(gra_intra[[4]],layout=layout.fruchterman.reingold)#edge.width=E(graj[[4]])$weight*10)
#adj = get.adjacency(gra_intra[[4]],attr='weight',sparse=FALSE)
if(FALSE){
plot(gra_intra[[4]],				#the graph to be plotted
     layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
     main='layer 4',	#specifies the title
     vertex.label.dist=0.5,			#puts the name labels slightly off the dots
     vertex.frame.color='blue', 		#the color of the border of the dots 
     vertex.label.color='black',		#the color of the name labels
     vertex.label.font=2,			#the font of the name labels
     vertex.label=NA, #V(graj[[4]])$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
     vertex.label.cex=1,			#specifies the size of the font of the labels. can also be made to vary
     vertex.size=5)
}

library(mvtnorm)
library(corrplot)
library(glmnet)
library(clusterGeneration)

k=10 # = Number of Candidate Variables
p=5 # = Number of Relevant Variables
N=500 # = Number of observations
betas=(-1)^(1:p) # = Values for beta
set.seed(12345) # = Seed for replication
sigma1=genPositiveDefMat(k,"unifcorrmat")$Sigma # = Sigma1 violates the irc
sigma2=sigma1 # = Sigma2 satisfies the irc
sigma2[(p+1):k,1:p]=0
sigma2[1:p,(p+1):k]=0

# = Verify the irrepresentable condition
irc1=sort(abs(sigma1[(p+1):k,1:p]%*%solve(sigma1[1:p,1:p])%*%sign(betas)))
irc2=sort(abs(sigma2[(p+1):k,1:p]%*%solve(sigma2[1:p,1:p])%*%sign(betas)))
c(max(irc1),max(irc2))

par(mfrow=c(1,2))
corrplot(cov2cor(sigma1))
corrplot(cov2cor(sigma2))
#type="upper", order="hclust", col=c("black", "white"),bg="lightblue")
data(iris)
classes <- iris$Species
variables <- iris[,1:4]
ccres <- corclust(variables, classes)
plot(ccres, mincor = 0.6)

library(igraph)
# prevent duplicated pairs
var.corelation <- cor(DATA[[1]][,sorted_index[1:100]+1])
check.corelation <- which(abs(var.corelation)>0.7, arr.ind=TRUE)

graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
groups.cor <- split(unique(as.vector(check.corelation)),clusters(graph.cor)$membership)
lapply(groups.cor,FUN=function(list.cor){rownames(var.corelation)[list.cor]})

library(igraph)
matt = DATA[[1]]
corr = cor(matt[,1:100])
colnames(corr) <- 1:40
check.corelation <- which(abs(corr)>0.7, arr.ind=TRUE)
graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
groups.cor <- split(unique(as.vector(check.corelation)),clusters(graph.cor)$membership)
lapply(groups.cor,FUN=function(list.cor){rownames(corr)[list.cor]})
out = corrplot(corr,order = "hclust", addrect = 8)
datas = dimnames(out)

corrMatOrder(corr, order = "hclust", hclust.method = "ward.D2")