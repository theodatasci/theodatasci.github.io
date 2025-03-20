library(cheddar)
library(igraph)
library(vegan)
library(sbm)
library(alluvial)
library(faux)
library(rmangal)
library(matlib)
library(calculus)
library(ggplot2)
#setwd("C:\\xxx")
source('functions_network.R')

####In FWebs, there is a large list of food webs called mg1
#mg1 <- create.fw.list(db="mg", ref=TRUE, spatial = TRUE)
load("mg1.Rdata")
sapply(1:length(mg1[[1]]),function(x) dim(as.matrix(mg1[[1]][[x]]))[1])

####Extract data from web 1,433
mat<-as.matrix(mg1[[1]][[433]])
plotMyMatrix(mat)
net.mat<-graph_from_adjacency_matrix(t(mat),mode="directed")

net<-sample_gnp(100,0.2,directed = FALSE)
plot(net)
net<-sample_gnp(100,0.2, directed =TRUE)
plot(net)

plot(net.mat,layout=layout_with_mds)
plot(net.mat,layout=layout_as_tree)
plot(net.mat,layout=layout_as_food_web5)

net<-sample_bipartite(50,50,"gnp",0.1)
plot(net,layout=layout_as_bipartite)

net<-sample_grg(100, 0.2)
get.adjacency(net)

E(net)$weight<-rpois(length(E(net)),2)
net[,]

degree(net)
degree(net,mode="in")
degree(net,mode="out")

####Compute the connectance of the empirical network
#answer
conn<-mean(as.matrix(net.mat[,]))

####Plot its degree distribution
#answer
hist(degree(net.mat),breaks=0:max(degree(net.mat)))
plot(degree_distribution(net.mat, cumulative = TRUE))

net<-sample_bipartite(4,4,"gnm",m=8)
net[,]

as_biadjacency_matrix(net)

m<-matrix(rbinom(9,1,0.5),nrow=3)
x<-rnorm(3)
m*x
m%*%x

eigen(m)

jacobian(c("x*y*z","x+y","z^2"), var = c("x", "y", "z"))

m<-matrix(rnorm(10^4),nrow=10^2)
plot(eigen(m)$values)

m<-cascade_matrix(10,20)
sum(m)/(dim(m)[1]*(dim(m)[1]-1))

m<-cascade_matrix(10,200)
sum(m)/(dim(m)[1]*(dim(m)[1]-1))

####Generate 100 virtual food webs based on the niche model with the actual connectance
#answer
niches<-lapply(1:100,function(x) niche_matrix(conn,dim(mat)[1]))
ms<-lapply(1:100,function(x) niches[[x]]$matrix)

niche<-niche_matrix(0.2,100)
m<-niche$matrix
sum(m)/(dim(m)[1]^2)

niche<-niche_matrix(0.2,200)
m<-niche$matrix
sum(m)/(dim(m)[1]^2)

net<-graph_from_adjacency_matrix(m,mode="directed")
i_index <- seq(from = 0, to = 1, by =0.1)
i_index <- head(i_index,-1)
prob_exp<-exponent.removal(net, i_index)
V(net)$name<-1:200
iterate(fw_to_attack=net, prob_exp, alpha1=50, iter=10, i_index, plot = TRUE)

net<-graph_from_adjacency_matrix(m)
hist(degree(net),breaks=0:max(degree(net)))
plot(degree_distribution(net, cumulative = TRUE))

net<-sample_degseq(degree(net),method = "vl")

mean(degree(net))
ks.test(degree(net),"pbinom",length(V(net)),mean(degree(net))/length(V(net)))

net<-sample_bipartite(50,50,"gnp",0.1)
sample.bip.config<-simulate(nullmodel(as_biadjacency_matrix(net),"curveball"),nsim=1000)
dim(sample.bip.config)

n<-1000
net<-sample_gnp(n,0.2, directed = TRUE)
sample.config.directed<-lapply(1:100,function(x) sample_degseq(degree(net,mode="out"), degree(net,mode="in"), method = "simple.no.multiple"))
length(sample.config.directed)

net<-sample_gnp(n,0.2, directed = FALSE)
sample.config.undirected<-lapply(1:100,function(x) sample_degseq(degree(net), method = "vl"))
length(sample.config.undirected)

net<-sample_gnp(100,0.2, directed = FALSE)
EB.mod<-cluster_edge_betweenness(net)
LE.mod<-cluster_leading_eigen(net)
ML.mod<-cluster_louvain(net)

plot(EB.mod,net,layout = layout_with_mds)
plot(LE.mod,net,layout = layout_with_mds)
plot(ML.mod,net,layout = layout_with_mds)

####Is the empirical FW less modular (or more modular) than the niche model food webs?
#answer
modul<-cluster_louvain(graph_from_adjacency_matrix(mat,mode="undirected"))
moduls<-sapply(1:100,function(x) cluster_louvain(graph_from_adjacency_matrix(ms[[x]],mode="undirected"))$modularity[2])
plot(density(moduls))
modul.ecdf<-ecdf(moduls)
1-modul.ecdf(modul$modularity[2])

p.val(modul$modularity[2],moduls,"larger","Null modularity distribution")

sbmnet <- sampleSimpleSBM(100, c(.5, .25, .25), list(mean = diag(.4, 3) + 0.05), model = 'bernoulli')
sbmnet$networkData
net.SBM <- estimateSimpleSBM(as.matrix(sbmnet$networkData))
plot(net.SBM, 'expected')
plot(net.SBM, 'data')

####is the empirical food web structured by modules or blocks?
#answer
m.SBM <- estimateSimpleSBM(mat)
make_alluvial_2(m.SBM$memberships,modul$membership,"Blocks","Modules")

####are blocks/modules related to trophic levels?
#answer
count_components(net.mat)
net.comp<-components(net.mat)
tl.1<-trophic_levels(largest_component(net.mat))
plot(largest_component(net.mat),layout = layout_as_food_web)
plot(tl.1~as.factor(m.SBM$indMemberships[which(net.comp$membership==1),]%*%(1:3)),xlab="SBM group",ylab="Trophic level")
plot(tl.1~as.factor(modul$membership[which(net.comp$membership==1)]),xlab="module",ylab="Trophic level")

####Spectral clustering
SC<-spectral_clustering(graph_from_adjacency_matrix(sbmnet$networkData),3)
plotMyMatrix(sbmnet$networkData,clustering=list("row"=SC,"col"=SC))
laplacian_spectral_gap(graph_from_adjacency_matrix(sbmnet$networkData))

optims<-sapply(4:10,function(n) laplacian_spectral_gap(graph_from_adjacency_matrix(sampleSimpleSBM(100, rep(1/n, n), list(mean = diag(.4, n) + 0.05), model = 'bernoulli')$networkData))$optim_n)
par(mfrow=c(1,1))
plot(4:10,optims)
