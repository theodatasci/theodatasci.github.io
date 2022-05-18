corrected_laplacian<-function(network){
	lap<-laplacian_matrix(network)
	lap+diag(degree(network,mode="in"))
}

trophic_levels<-function(network) { #based on MacKay et al. 2020
	lap<-corrected_laplacian(network)
	imbalance<-degree(network,mode="in")-degree(network,mode="out")
	inv(as.matrix(lap)) %*% imbalance
}

cascade_matrix<-function(cc,nspecies){ #from Cohen-Newman-Briand's papers
	upper.tri(matrix(1,nrow=nspecies,ncol=nspecies), diag = FALSE)*matrix(rbinom(nspecies*nspecies,1,cc/nspecies),nrow=nspecies)
}

niche_matrix<-function(connectance,nspecies){ #Williams-Martinez model
	n_i<-runif(nspecies)
	r_i<-rbeta(nspecies,1,(1/(2*connectance))-1)
	r_i<-r_i*n_i
	c_i<-sapply(1:nspecies,function(x) runif(1,min=r_i[x]/2,max=n_i[x]))
	pred_function<-function(z_1,z_2){
		if((n_i[z_1]<=c_i[z_2]+0.5*r_i[z_2])&(n_i[z_1]>=c_i[z_2]-0.5*r_i[z_2])) {
			1
		}
		else{
			0
		}
	}
	mat<-sapply(1:nspecies,function(x) sapply(1:nspecies,function(y) pred_function(y,x)))
	list("matrix"=mat,"n"=n_i,"r"=r_i,"c"=c_i)
}


to_upper_triangular<-function(mat){
	upper.tri(mat, diag = FALSE)*mat
}

make_alluvial_2<-function(classif1,classif2,name1,name2){
	A <- as.data.frame(table(classif1,classif2))
	colnames(A) = c(name1,name2,"Freq")
	w   <- which(A$Freq != 0)
	A <- A[w,]
	alluvial(A[,c(1,2)],freq = A$Freq)
}

FW_interaction_from_predation<-function(mat,rho){
	n<-dim(mat)[1]
	fill<-rnorm_multi(n*(n-1)/2,2,r=rho)
	ut<-matrix(0,n,n)
	ut[lower.tri(ut,diag = FALSE)]<-fill[,1]
	ut<-mat*t(ut)
	lt<-matrix(0,n,n)
	lt[lower.tri(lt, diag = FALSE)]<-fill[,2]
	lt<-(t(mat))*lt
	ut+lt
}
