setwd("C:\\Massol\\Enseignement\\Théorie CESAB")

library(igraph)
library(matlib)
source('functions_network.R')
n<-100

#generate random unidrected unipartite network
net<-erdos.renyi.game(n,0.2, "gnp")
plot(net)

#generate random directed unipartite network
net<-erdos.renyi.game(n,0.2, "gnp",directed =TRUE)
plot(net)
V(net)$y<-trophic_levels(net)
V(net)$x<-layout_with_graphopt(net)[,1]
plot(net,layout=layout_nicely)

#generate random bipartite network
net<-sample_bipartite(n/2,n/2,"gnp",0.1)
plot(net,layout=layout_as_bipartite)

net<-sample_grg(100, 0.2)
get.adjacency(net)

E(net)$weight<-rpois(length(E(net)),2)
net[,]

#degrees
degree(net)
net<-erdos.renyi.game(n,0.2, "gnp",directed =TRUE)
degree(net,mode="in")
degree(net,mode="out")

#adjacency vs. incidence matrices
net<-sample_bipartite(4,4,"gnm",m=8)
net[,]
as_incidence_matrix(net)

#matrix multiplication
m<-matrix(rbinom(9,1,0.5),nrow=3)
x<-rnorm(3)
m*x
m%*%x

#matrix spectra
eigen(m)

#computing Jacobians
library(calculus)
jacobian(c("x*y*z","x+y","z^2"), var = c("x", "y", "z"))

#May's result
m<-matrix(rnorm(10^4),nrow=10^2)
plot(eigen(m)$values)
m<-matrix(rnorm(10^6),nrow=10^3)
plot(eigen(m)$values)

#cascade model
m<-cascade_matrix(10,20)
sum(m)/(dim(m)[1]*(dim(m)[1]-1))
m<-cascade_matrix(10,200)
sum(m)/(dim(m)[1]*(dim(m)[1]-1))

#niche model
niche<-niche_matrix(0.2,100)
m<-niche$matrix
sum(m)/(dim(m)[1]^2)

niche<-niche_matrix(0.2,200)
m<-niche$matrix
sum(m)/(dim(m)[1]^2)

#robustness analysis see https://github.com/FMestre1/fw_package
#library(devtools) 
#install_github("FMestre1/fw_package")
library(FWebs)
net<-graph_from_adjacency_matrix(m,mode="directed")
i_index <- seq(from = 0, to = 1, by =0.1)
i_index <- head(i_index,-1)
prob_exp<-exponent.removal(net, i_index)
V(net)$name<-1:200
iterate(fw_to_attack=net, prob_exp, alpha1=50, iter=10, i_index, plot = TRUE)

#degree distribution
net<-graph_from_adjacency_matrix(m)
hist(degree(net),breaks=0:max(degree(net)))
plot(degree_distribution(net, cumulative = TRUE))

#generate unipartite network from degree sequence
net<-sample_degseq(degree(net),method = "vl")

#compare degree distributions
mean(degree(net))
ks.test(degree(net),"pbinom",length(V(net)),mean(degree(net))/length(V(net)))

#randomization of bipartite network (Strona's curveball algorithm)
library(vegan)
net<-sample_bipartite(50,50,"gnp",0.1)
sample.bip.config<-simulate(nullmodel(as_incidence_matrix(net),"curveball"),nsim=1000)
dim(sample.bip.config)

#randomization of directed networks
net<-erdos.renyi.game(n,0.2, "gnp",directed =TRUE)
sample.config.directed<-lapply(1:100,function(x) sample_degseq(degree(net,mode="out"), degree(net,mode="in"), method = "simple.no.multiple"))
length(sample.config.directed)

#randomization of unipartite networks
net<-erdos.renyi.game(n,0.2, "gnp")
sample.config.undirected<-lapply(1:100,function(x) sample_degseq(degree(net), method = "vl"))
length(sample.config.undirected)

#modularity
net<-erdos.renyi.game(100,0.2, "gnp",directed =FALSE)
EB.mod<-cluster_edge_betweenness(net)
LE.mod<-cluster_leading_eigen(net)
ML.mod<-cluster_louvain(net)

plot(EB.mod,net,layout = layout_with_mds)
plot(LE.mod,net,layout = layout_with_mds)
plot(ML.mod,net,layout = layout_with_mds)

#block model
library(sbm)
sbmnet <- sampleSimpleSBM(100, c(.5, .25, .25), list(mean = diag(.4, 3) + 0.05), model = 'bernoulli')
sbmnet$networkData
net.SBM <- estimateSimpleSBM(as.matrix(sbmnet$networkData))
plot(net.SBM, 'expected')
plot(net.SBM, 'data')
