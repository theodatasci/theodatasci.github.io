library(igraph)
library(vegan)
library(sbm)
library(alluvial)
library(faux)
library(FWebs)
source('functions_network.R')

####In FWebs, there is a large list of food webs called mg1
sapply(1:length(mg1[[1]]),function(x) dim(as.matrix(mg1[[1]][[x]]))[1])

####Extract data from web 1,223
mat<-as.matrix(mg1[[1]][[223]])
plotMyMatrix(mat)
net<-graph_from_adjacency_matrix(mat,mode="directed")


####Q0: Compute the connectance of the network

####Q1: Plot the degree distribution

####Q2: is this food web structured by modules or blocks?

####Q3: are blocks/modules related to trophic levels?

####Q4: Generate 100 virtual food webs based on the niche model with the actual connectance

####Q5: Is the actual FW less modular (or more modular) than the niche model food webs?


