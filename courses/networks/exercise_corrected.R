library(igraph)
library(vegan)
library(sbm)
library(alluvial)
library(faux)
library(FWebs)
library(matlib)
source('functions_network.R')

####In FWebs, there is a large list of food webs called mg1
sapply(1:length(mg1[[1]]),function(x) dim(as.matrix(mg1[[1]][[x]]))[1])

####Extract data from web 1,223
mat<-as.matrix(mg1[[1]][[223]])
plotMyMatrix(mat)
net<-graph_from_adjacency_matrix(mat,mode="directed")


####Q0: Compute the connectance of the network
conn<-mean(as.matrix(net[,]))

####Q1: Plot the degree distribution
hist(degree(net),breaks=0:max(degree(net)))
plot(degree_distribution(net, cumulative = TRUE))


####Q2: is this food web structured by modules or blocks?
m.SBM <- estimateSimpleSBM(mat)
modul<-cluster_louvain(graph_from_adjacency_matrix(mat,mode="undirected"))
make_alluvial_2(m.SBM$memberships,modul$membership,"Blocks","Modules")


####Q3: are blocks/modules related to trophic levels?
count_components(net)
net.comp<-components(net)
tl.1<-trophic_levels(induced_subgraph(net, which(net.comp$membership==1)))

plot(tl.1~as.factor(m.SBM$indMemberships[which(net.comp$membership==1),]%*%(1:3)),xlab="SBM group",ylab="Trophic level")
plot(tl.1~as.factor(modul$membership[which(net.comp$membership==1)]),xlab="module",ylab="Trophic level")


####Q4: Generate 100 virtual food webs based on the niche model with the actual connectance
niches<-lapply(1:100,function(x) niche_matrix(conn,dim(mat)[1]))
ms<-lapply(1:100,function(x) niches[[x]]$matrix)



####Q5: Is the actual FW less modular (or more modular) than the niche model food webs?
moduls<-lapply(1:100,function(x) cluster_louvain(graph_from_adjacency_matrix(ms[[x]],mode="undirected"))$modularity)
plot(density(unlist(moduls)))

modul.ecdf<-ecdf(unlist(moduls))
modul.ecdf(modul$modularity[2])



















#TRASH
dim(as.matrix(mg1[[1]][[223]]))
plot(m.SBM, 'expected')
plot(m.SBM, 'data')


sapply(1:length(mg1), function(x) sapply(1:length(mg1[[x]]), function(y) dim(as.matrix(mg1[[x]][[y]]))[1]))




plotMyMatrix(ms[[1]][order(niches[[1]]$n),order(niches[[1]]$n)])
modul$modularity[2]

net.niche<-graph_from_adjacency_matrix(m,mode="directed")
count_components(net.niche)


plotMyMatrix(FW_interaction_from_predation(mat,-0.6))
plotMyMatrix(FW_interaction_from_predation(matrix(rbinom(100,1,0.5),ncol=100,nrow=100),-0.6))








niche<-niche_matrix(0.2,100)
m<-niche$matrix
signs<-m-t(m)
sym<-m+t(m)
diag(sym)<-rep(0,100)
jac<-signs*rnorm(100*100,0,1/100)
jac<-jac-diag(rep(1,100))
sp<-eigen(jac)
plot(sp$value)

m.SBM <- estimateSimpleSBM(m)
plot(m.SBM, 'data')
m2<-m[order(niche$n),order(niche$n)]
plotMyMatrix(m2)

net<-graph_from_adjacency_matrix(m,mode="directed")
tl<-trophic_levels(net)

V(net)$y<-niche$n
V(net)$x<-niche$c
plot(net,layout=layout_nicely)

plot(tl~niche$n)
summary(lm(tl~niche$n))$adj.r.squared

#configs.directed<-lapply(1:100,function(x) sample_degseq(degree(net,mode="out"), degree(net,mode="in"), method = "simple.no.multiple"))
configs.undirected<-lapply(1:100,function(x) sample_degseq(sapply(1:100,function(y) sum(sym[y,])), method = "vl"))
configs.permutations<-lapply(1:100, function(x) sample(1:100))
configs.undirected.UT<-lapply(1:100,function(x) to_upper_triangular(as.matrix(configs.undirected[[x]][configs.permutations[[x]],configs.permutations[[x]]])))
configs.adj.r.squared<-lapply(1:100,function(x) summary(lm(trophic_levels(graph_from_adjacency_matrix(configs.undirected.UT[[x]],mode="directed"))~(niche$n)[configs.permutations[[x]]]))$adj.r.squared)

plot(density(unlist(configs.adj.r.squared)))
ecdf.adj.r.squared<-ecdf(unlist(configs.adj.r.squared))
1-ecdf.adj.r.squared(summary(lm(tl~niche$n))$adj.r.squared)



configs.undirected.UT.signs<-lapply(1:100, function(x) matrix((2*rbinom(100*100,1,0.5)-1),nrow=100)*configs.undirected.UT[[x]])
configs.undirected.signs<-lapply(1:100, function(x) configs.undirected.UT.signs[[x]] - t(configs.undirected.UT.signs[[x]]))


i_index <- seq(from = 0, to = 1, by =0.01)
i_index <- head(i_index,-1)
prob_exp<-exponent.removal(net, i_index)
V(net)$name<-1:100
iterate(fw_to_attack=net, prob_exp, alpha1=50, iter=10, i_index, plot = TRUE)




