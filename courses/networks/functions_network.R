corrected_laplacian<-function(network){
	lap<-laplacian_matrix(network)
	lap+diag(degree(network,mode="in"))
}

trophic_levels<-function(network) { #based on MacKay et al. 2020
	lap<-corrected_laplacian(network)
	imbalance<-degree(network,mode="in")-degree(network,mode="out")
	inv(as.matrix(lap)) %*% imbalance
}

as_Community<-function(network,net_name="."){#takes a directed network (with species names) and make it a Community object for cheddar
	names_net<-V(network)$name
	nodes<-data.frame("node"=names_net,row.names= names_net)
	properties=list("title"=net_name)
	tlinks<-as_long_data_frame(network)
	trophic.links<-as.matrix(tlinks[,c("from_name","to_name")])
	colnames(trophic.links)<-c("resource", "consumer")
	Community(nodes, properties, trophic.links)
}
#example: ShortestTrophicLevel(as_Community(network,"WTF"))

layout_as_food_web<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-TrophicLevels(as_Community(network))[,5]-1
	lay
}

layout_as_food_web2<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-trophic_levels(network)
	lay
}

layout_as_food_web3<-function(network){#adapted from Jon Borelli's code https://assemblingnetwork.wordpress.com/2013/04/03/more-food-web-plotting-with-r/
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,1]<-layout_with_graphopt(network)[,1]
	lay[,2]<-alt_TL(network)[,3]
	lay
}

layout_as_food_web4<-function(network){
	graphlayouts::layout_with_constrained_stress(network,coord = alt_TL(network)[,3],fixdim = "y")
}

layout_as_food_web5<-function(network){
	l<-length(V(network))
	lay<-matrix(nrow=l,ncol=2) 
	lay[,2]<-alt_TL(network)[,3]
	n<-table(lay[,2])
	vals<-sort(unique(lay[,2]))
	for(i in 1:length(n)){
		locations<-which(lay[,2]==vals[i])
		lay[locations,1]<-(1:n[i])/n[i]
	}
	lay
}

alt_TL<-function(network){ #yet another implementation of trophic levels
	undir_net<-as.undirected(network)
	basals<-which(degree(network,mode="in")==0)
	dist_mat<-t(distances(network,v=which(degree(network,mode="in")==0),mode="out"))
	s<-dim(dist_mat)[1]
	shortest_chain<-sapply(1:s,function(x) min_without_inf(dist_mat[x,]))+1
	longest_chain<-sapply(1:s,function(x) max_without_inf(dist_mat[x,]))+1
	average_chain<-sapply(1:s,function(x) mean_without_inf(dist_mat[x,]))+1
	if(!is.null(V(network)$name)){
		res<-data.frame("species"=V(network)$name,"shortest"=shortest_chain,"longest"=longest_chain,"average" = average_chain)
	}
	else{
		res<-data.frame("shortest"=shortest_chain,"longest"=longest_chain,"average" = average_chain)
	}
	res
}

min_without_inf<-function(vec){
	min(vec[!is.infinite(vec)])
}

max_without_inf<-function(vec){
	max(vec[!is.infinite(vec)])
}

mean_without_inf<-function(vec){
	mean(vec[!is.infinite(vec)])
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

spectral_clustering <- function(graph, nb_cluster, normalized = TRUE) {#from J. Chiquet's git page https://jchiquet.github.io/MAP566/docs/mixture-models/map566-lecture-graph-clustering-part1.html
  
  ## Compute Laplacian matrix
  L <- igraph::laplacian_matrix(graph, normalized = normalized) 
  ## Generates indices of last (smallest) K vectors
  selected <- rev(1:ncol(L))[1:nb_cluster] 
  ## Extract n normalized eigen-vectors
  U <- eigen(L)$vectors[, selected, drop = FALSE]  # spectral decomposition
  U <- sweep(U, 1, sqrt(rowSums(U^2)), '/')    
  ## Perform k-means
  res <- kmeans(U, nb_cluster, nstart = 40)$cl
  
  res
}

laplacian_spectral_gap<- function(graph){
	L <- igraph::laplacian_matrix(graph, normalized = TRUE)
	comps<-igraph::count_components(graph)
	lambdas <- sort(eigen(L)$values)
	l<-length(lambdas)
	s_gaps<-lambdas[-1]-lambdas[-l]
	s_util<-s_gaps[-(1:comps)]
	s_util<-s_util[1:round(l/2)]
	opt_n<-which.max(s_util)+comps
	
	par(mfrow=c(2,1))
	plot(lambdas,xlab="",ylab="lambda",type="l")
	plot(s_gaps,xlab="",ylab="spectral gap",type="l")
	
	list("spectral_gaps"=s_gaps,"optim_n"=opt_n)
}

p.val<-function(test_val,test_collection,method="larger",label=""){#compute a p-value for a test based on a collection of values simulated using a null model
	test_collection<-c(test_collection,test_val)
	n<-length(test_collection)
	cumul<-ecdf(test_collection)
	mmin<-min(test_collection)
	mmax<-max(test_collection)
	plot(density(test_collection),xlim=c(mmin,mmax),main=label,xlab="")
	abline(v=test_val,col="red")
	if(method=="lower"){
		max(cumul(test_val),1/n)
	}
	else if(method=="two-sided"){
		max(2*min(cumul(test_val),1-cumul(test_val)),1/n)
	}
	else{
		max(1-cumul(test_val),1/n)
	}
	
}

exponent.removal<-function (fw, i_index){
    d1 <- as.data.frame(degree(fw))
    kmax <- max(d1$`degree(fw)`)
    kmin <- min(d1$`degree(fw)`)
    k_vector <- sort(unique(d1$`degree(fw)`))
    colnames_Pe <- paste0("Pe_", i_index)
    Pe_results <- as.data.frame(matrix(ncol = length(colnames_Pe), 
        nrow = nrow(d1)))
    Pe_results <- cbind(d1, Pe_results)
    colnames(Pe_results)[-1] <- colnames_Pe
    for (j in 1:nrow(d1)) {
        current_species_k <- d1[j, ]
        Nk <- sum(d1$`degree(fw)` == current_species_k)
        for (i in 1:length(i_index)) {
            current_i_value <- i_index[i]
            sum_denominator <- as.data.frame(matrix(nrow = length(k_vector)))
            sum_denominator <- cbind(k_vector, sum_denominator)
            for (f in 1:nrow(sum_denominator)) {
                Ni <- sum(d1$`degree(fw)` == sum_denominator[f, 
                  1])
                sum_denominator[f, 2] <- ((1 - current_i_value)^(kmax - 
                  sum_denominator[f, 1])) * Ni
            }
            sum_denominator <- sum(sum_denominator$V1)
            Pe <- (((1 - current_i_value)^(kmax - current_species_k)) * 
                Nk)/sum_denominator
            Pe_results[j, i + 1] <- Pe
        }
    }
    return(Pe_results[, -1])
}

iterate<-function (fw_to_attack, probs_of_fw, alpha1, iter, i_index, plot = FALSE, export_plot = FALSE, plot_name = NULL){
    result_iterate <- data.frame(matrix(nrow = ncol(probs_of_fw)))
    for (i in 1:iter) {
        r1 <- robustness(fw_to_attack, probs_of_fw, alpha1 = 50)
        R1 <- r1$Ralpha
        result_iterate <- cbind(result_iterate, R1)
        message(paste0("Iteration ", i))
    }
    result_iterate <- result_iterate[, -1]
    meanR <- apply(result_iterate, 1, FUN = mean)
    sdR <- apply(result_iterate, 1, FUN = sd)
    output <- as.data.frame(cbind(i_index, meanR, sdR))
    output.se <- output$sdR/sqrt(nrow(output))
    margin.error <- qnorm(0.975) * output.se
    lower.bound <- output$meanR - margin.error
    upper.bound <- output$meanR + margin.error
    output <- data.frame(output, lower.bound, upper.bound)
    if (any(output$lower.bound < 0)) 
        output[output$lower.bound < 0, ]$lower.bound <- 0
    if (plot == TRUE) {
        print(ggplot(output, aes(x = i_index, y = meanR), xlab = "label") + 
            xlab("Intentionality (I)") + ylab(paste0("R", alpha1)) + 
            ylim(0, (alpha1/100) + 0.1) + geom_line(color = "steelblue4", 
            lwd = 1) + geom_ribbon(alpha = 0.5, aes(ymin = lower.bound, 
            ymax = upper.bound), fill = "steelblue2", color = "steelblue2"))
    }
    if (export_plot == TRUE) {
        png(paste0(plot_name, ".png"), width = 500, height = 400)
        print(ggplot(output, aes(x = i_index, y = meanR), xlab = "label") + 
            xlab("Intentionality (I)") + ylab(paste0("R", alpha1)) + 
            ylim(0, (alpha1/100) + 0.1) + geom_line(color = "steelblue4", 
            lwd = 1) + geom_ribbon(alpha = 0.5, aes(ymin = lower.bound, 
            ymax = upper.bound), fill = "steelblue2", color = "steelblue2"))
        dev.off()
    }
    return(output)
}

dd.fw<-function (list1, log = TRUE, cumulative = TRUE) {
    if (class(list1[[1]]) == "list") {
        list2 <- list1[[1]]
    }
    else list2 <- list1
    df_final <- data.frame()
    for (i in 1:length(list2)) {
        m1 <- list2[[i]]
        m1 <- as.matrix(m1)
        g1 <- igraph::graph_from_adjacency_matrix(m1, weighted = NULL, 
            mode = "undirected")
        g2 <- igraph::degree_distribution(g1, cumulative)
        d <- igraph::degree(g1, mode = "all")
        degree1 <- 1:max(d)
        probability <- g2[-1]
        nonzero.position <- which(probability != 0)
        probability <- (probability[nonzero.position])
        degree1 <- (degree1[nonzero.position])
        iter <- rep(paste0("iter", i), length(degree1))
        colour1 <- rep(randomColor(), length(iter))
        df0 <- data.frame(iter, degree1, probability, colour1)
        df_final <- rbind(df_final, df0)
    }
    if (log == TRUE) {
        print(ggplot2::ggplot(df_final, aes(x = degree1, y = probability, 
            group = iter)) + geom_line(aes(col = factor(iter))) + 
            theme(legend.position = "none") + labs(title = "Degree distribution (log-log)", 
            x = "degree (log)", y = "Probability (log)") + theme(plot.title = element_text(hjust = 0.5)) + 
            scale_y_log10() + scale_x_log10() + annotation_logticks())
    }
    if (log == FALSE) {
        print(ggplot2::ggplot(df_final, aes(x = degree1, y = probability, 
            group = iter)) + geom_line(aes(col = factor(iter))) + 
            theme(legend.position = "none") + labs(title = "Degree distribution", 
            x = "degree", y = "Probability") + theme(plot.title = element_text(hjust = 0.5)))
    }
    return(df_final[, -4])
}

create.fw.list<-function (db, folder = NULL, ecosyst = FALSE, ref = FALSE, spatial = FALSE, code = FALSE, mangal_types = NULL){
    model.dissemination_allow <- model.whole_food_web <- NULL
    fwlist <- list()
    if (!db %in% c("eb", "wl", "gw", "mg")) 
        stop("Argument 'db' must take one of the following values:\n\n                                          'wl' - Web of Life\n                                          'mg' - mangal\n                                          'gw' - globalweb\n                                          'eb' - ecobase")
    if (!db %in% c("wl", "gw") & !is.null(folder)) 
        stop("Argument 'folder'can only be used if 'db'= 'wl' or 'gw'!")
    if (!db %in% c("mg") & !is.null(mangal_types)) 
        stop("Argument 'type'can only be used if 'db'= 'mg'!")
    if (!db %in% c("gw", "eb") & ecosyst == TRUE) 
        stop("Argument 'ecosyst'can only be used if 'db'= 'eb' or 'gw'!")
    if (!db %in% c("wl", "mg", "eb") & spatial == TRUE) 
        stop("Argument 'spatial'can only be used if 'db'= 'eb', 'mg' 'wl'!")
    if (!db %in% c("wl", "mg", "gw") & code == TRUE) 
        stop("Argument 'code'can only be used if 'db'= 'wl', 'mg'', 'gw'!")
    if ("mg" %in% db & is.null(mangal_types)) 
        message("No value defined for the 'mangal_types' argument! \n Will assume types 'predation' and 'herbivory'.")
    if (!"mg" %in% db & !is.null(mangal_types)) 
        stop("Argument 'mangal_types'can only be used if 'db'= 'mg'!")
    if (db == "gw") {
        message("####################### GLOBALWEB DATABASE #######################\n\n")
        message("Fetching info from the provided folder!")
        files_gw <- list.files(path = folder, pattern = "WEB")
        ngw <- length(files_gw)
        message(paste0("There are ", ngw, " food web files in the folder!"))
        message("You should have downloaded the file 'Current Food Web List' from the GlobalWeb website\n             \n and converted it to csv.")
        if (ref == TRUE) 
            reflist_gw <- c()
        names_gw <- c()
        for (i in 1:ngw) {
            message(paste0("Fetching food web ", i, " in ", ngw, 
                "!"))
            dfgw <- read.csv(paste0(folder, "/", files_gw[i]), 
                header = FALSE)
            dfgw <- dfgw[, colSums(is.na(dfgw)) <= 1]
            names_gw[i] <- as.character(dfgw[2, 1])
            if (ref == TRUE) 
                reflist_gw[i] <- as.character(dfgw[1, 1])
            names_gw_c <- c()
            n1 <- ncol(dfgw) - 1
            for (j in 1:n1) {
                names_gw_c[j] <- as.character(dfgw[2, j + 1])
            }
            names_gw_r <- c()
            n2 <- nrow(dfgw) - 2
            for (j in 1:n2) {
                names_gw_r[j] <- as.character(dfgw[j + 2, 1])
            }
            dfgw <- dfgw[-c(1, 2), -1]
            dfgw[dfgw == ""] <- NA
            dfgw <- na.omit(dfgw)
            if (i == 281) {
                names_gw_r <- names_gw_r[-c(36, 37)]
            }
            names_gw_c <- names_gw_c[names_gw_c != ""]
            names_gw_r <- names_gw_r[names_gw_r != ""]
            names_gw_c <- paste0("sp_", as.character(1:length(names_gw_c)), 
                "_", names_gw_c)
            names_gw_r <- paste0("sp_", as.character(1:length(names_gw_r)), 
                "_", names_gw_r)
            colnames(dfgw) <- names_gw_c
            rownames(dfgw) <- names_gw_r
            fwlist[[i]] <- dfgw
        }
        names(fwlist) <- names_gw
        if (ref == TRUE) {
            references <- as.data.frame(matrix(ncol = 4))
            names(references) <- c("FW code", "first_author", 
                "year", "full_ref")
            files_gw <- list.files(folder, pattern = "WEB")
            message("Fetching references from the dataset files!")
            for (w in 1:ngw) {
                dfgw <- read.csv(paste0(folder, "/", files_gw[w]), 
                  header = FALSE)
                dfgw <- dfgw[, colSums(is.na(dfgw)) <= 1]
                full_ref1 <- as.character(dfgw[1, 1])
                references[w, 4] <- full_ref1
                references[w, 1] <- files_gw[w]
                references[w, 2] <- str_sub(word(full_ref1, start = 1), 
                  1, str_length(word(full_ref1, start = 1)) - 
                    1)
                references[w, 3] <- regmatches(x = full_ref1, 
                  gregexpr("[0-9]+", text = full_ref1))[[1]][1]
            }
        }
        if (ecosyst == TRUE) {
            message("Searching for 'gw_list.csv' file...")
            if (!file.exists(paste0(folder, "/gw_list.csv"))) 
                stop("\nDownload the file 'Current Food Web List' from the website\n                                                             \nand convert to a csv named 'gw_list.csv' please!")
            gw_eco <- read.csv(paste0(folder, "/", "gw_list.csv"), 
                header = TRUE, sep = ";")
            filn <- paste0("WEB", as.character(gw_eco[, 1]), 
                ".csv")
            gw_eco2 <- gw_eco[, 1:3]
            gw_eco2[, 1] <- filn
            names(gw_eco2)[1] <- "FW"
            filn <- as.data.frame(cbind(filn, filn))
            names(filn) <- c("filn1", "filn2")
            ecosystem <- merge(x = filn, y = gw_eco2, by.x = "filn2", 
                by.y = "FW")
            ecosystem <- ecosystem[, c(2, 3, 4)]
            names(ecosystem)[1] <- "Food web"
        }
    }
    if (db == "wl") {
        message("####################### WEB OF LIFE DATABASE #######################\n\n")
        files_wl <- list.files(path = folder, pattern = ".csv")
        files_wl <- files_wl[files_wl != "references.csv"]
        nwl <- length(files_wl)
        message(paste0("There are ", nwl, " food web files in the folder!"))
        if (file.exists(paste0(folder, "/references.csv"))) {
            table_wl <- read.csv(paste0(folder, "/references.csv"), 
                header = TRUE)
        }
        else {
            stop("There is no 'references.csv' file on the folder, as provided by the website!")
        }
        names_wl <- as.character(table_wl[, 8])
        for (i in 1:nwl) {
            message(paste0("Fetching food web ", i, " in ", nwl, 
                "!"))
            dfwl <- read.csv(paste0(folder, "/", files_wl[i]), 
                header = TRUE, row.names = 1)
            dfwl[is.na(dfwl)] <- 0
            fwlist[[i]] <- dfwl
        }
        names(fwlist) <- names_wl
        if (ref == TRUE) {
            references <- as.data.frame(matrix(ncol = 4))
            names(references) <- c("FW code", "first_author", 
                "year", "full_ref")
            message("Fetching references from the 'references.csv' file!")
            message("Checking the presence of the 'references.csv' file...")
            if (!file.exists(paste0(folder, "/references.csv")) == 
                TRUE) 
                stop("Can't retrieve reference details... \n File not present!")
            ref_file <- read.csv(paste0(folder, "/references.csv"), 
                header = TRUE)
            for (w in 1:nwl) {
                full_ref1 <- as.character(ref_file[w, 7])
                references[w, 4] <- full_ref1
                references[w, 1] <- as.character(ref_file[w, 
                  1])
                references[w, 2] <- stringr::str_sub(stringr::word(full_ref1, 
                  start = 1), 1, stringr::str_length(stringr::word(full_ref1, 
                  start = 1)) - 1)
                references[w, 3] <- regmatches(x = full_ref1, 
                  gregexpr("[0-9]+", text = full_ref1))[[1]][1]
            }
        }
        if (spatial == TRUE) {
            message("Fetching the spatial information from the 'references.csv' file!")
            message("Checking the presence of the 'references.csv' file...")
            if (!file.exists(paste0(folder, "/references.csv")) == 
                TRUE) 
                stop("Can't retrieve spatial info... \n File not present!")
            ref_file <- read.csv(paste0(folder, "/references.csv"), 
                header = TRUE)
            spatial1 <- ref_file[, c(1, 9, 10)]
        }
    }
    if (db == "eb") {
        message("####################### ECOBASE DATABASE #######################\n\n")
        message("Fetching info from the EcoBase website!")
        suppressWarnings({
            suppressMessages({
                h = basicTextGatherer()
                curlPerform(url = "http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client_3.php", 
                  writefunction = h$update)
                data1 <- xmlTreeParse(h$value(), useInternalNodes = TRUE)
                liste_mod <- ldply(xmlToList(data1), data.frame)
            })
            l2 <- subset(liste_mod, model.dissemination_allow == 
                "true")
            message("Sellected only those to which model dissemination is allowed!")
            l3 <- subset(l2, model.whole_food_web == "true")
            message("Sellected only those to which the whole food web is available!")
            model.name <- as.character(l3$model.model_name)
            input_list <- list()
            id <- as.numeric(as.character(l3$model.model_number))
            for (i in 1:nrow(l3)) {
                message(paste0("Fetching information on food web ", 
                  i, " of ", nrow(l3)))
                suppressMessages({
                  h = basicTextGatherer()
                  mymodel <- id[i]
                  curlPerform(url = paste("http://sirs.agrocampus-ouest.fr/EcoBase/php/webser/soap-client.php?no_model=", 
                    mymodel, sep = ""), writefunction = h$update, 
                    verbose = TRUE)
                  data2 <- xmlTreeParse(h$value(), useInternalNodes = TRUE)
                  input1 <- xpathSApply(data2, "//group", function(x) xmlToList(x))
                })
                names_input <- as.character(input1[1, ])
                input1 <- as.data.frame(input1)
                colnames(input1) <- names_input
                input1 <- input1[-1, ]
                input_list[[i]] <- input1
            }
            mnames <- names(input_list)
            for (i in 1:length(input_list)) {
                m2 <- input_list[[i]]
                nnodes <- length(m2)
                node_names <- names(m2)
                int_matrix <- as.data.frame(matrix(ncol = nnodes, 
                  nrow = nnodes))
                for (j in 1:length(m2)) {
                  node1 <- m2[[j]]
                  node_id <- as.numeric(node1$group_seq)
                  node_name <- node_names[j]
                  colnames(int_matrix)[node_id] <- node_name
                  rownames(int_matrix)[node_id] <- node_name
                  diet_node1 <- node1$diet_descr
                  nr_food_items <- length(diet_node1)
                  for (a in 1:nr_food_items) {
                    item1 <- diet_node1[[a]]
                    id_item1 <- as.numeric(item1$prey_seq)
                    proportion_item1 <- as.numeric(item1$proportion)
                    detritus_item1 <- as.numeric(item1$detritus_fate)
                    int_matrix[id_item1, node_id] <- proportion_item1
                  }
                }
                int_matrix[is.na(int_matrix)] <- 0
                fwlist[[i]] <- int_matrix
            }
            names(fwlist) <- model.name
        })
        if (ref == TRUE) {
            references <- as.data.frame(matrix(ncol = 4))
            names(references) <- c("FW code", "first_author", 
                "year", "full_ref")
            message("Fetching the references information!")
            for (w in 1:nrow(l3)) {
                full_ref1 <- as.character(l3$model.reference)[w]
                references[w, 4] <- full_ref1
                references[w, 1] <- as.numeric(as.character(l3$model.model_number[w]))
                references[w, 2] <- as.character(l3$model.author[w])
                references[w, 3] <- regmatches(x = full_ref1, 
                  gregexpr("[0-9]+", text = full_ref1))[[1]][1]
            }
        }
        if (ecosyst == TRUE) {
            ecosystem <- data.frame(l3$model.model_number, l3$model.country, 
                l3$model.ecosystem_type)
            names(ecosystem) <- c("Food web", "Location", "Ecosystem")
        }
        if (spatial == TRUE) {
            message("Fetching spatial information from the EcoBase website...")
            if (!file.exists("ecobase_areas.shp")) {
                stop("If you need the spatial information on each dataset you have to:\n\n             1. Download the kml file from http://sirs.agrocampus-ouest.fr/EcoBase/php/protect/extract_kml.php;\n\n             (file name is 'location.kml')\n\n             2. Convert it to a shapefile in any GIS;\n\n             3. Name it 'ecobase_areas.shp';\n\n             4. Place it in the working directory;\n\n             ... I know, it is not ideal!...\n             ")
            }
            else EcoBase_shape <- sf::st_read("ecobase_areas.shp")
            ebd <- EcoBase_shape$Name
            nmr <- list()
            for (i in 1:length(ebd)) {
                nr <- strsplit(as.character(ebd[i]), "--::")[[1]][1]
                nr <- as.numeric(str_extract_all(nr, "\\d+")[[1]])
                nmr[[i]] <- nr
            }
            nmr2 <- c()
            for (i in 1:length(nmr)) {
                a <- nmr[[i]]
                b <- length(a)
                c1 <- rep(i, b)
                nmr2 <- c(nmr2, c1)
            }
            nmr <- unlist(nmr)
            table1 <- as.data.frame(cbind(nmr2, nmr))
            colnames(table1) <- c("row_n", "id")
            lines_n <- c()
            for (i in 1:nrow(liste_mod)) {
                id <- as.numeric(as.character(liste_mod$model.model_number[i]))
                lines_n[i] <- as.numeric(table1[table1$id == 
                  id, ][1])
            }
            ecobase_poly2 <- list()
            for (i in 1:length(lines_n)) {
                ecobase_poly2[i] <- st_geometry(EcoBase_shape)[lines_n[i]]
            }
            for (i in 1:length(ecobase_poly2)) {
                if (is.na(lines_n[i])) {
                  z1 <- as.numeric(unlist(regmatches(liste_mod$model.geographic_extent[[i]], 
                    gregexpr("[[:digit:]]+\\.*[[:digit:]]*", 
                      liste_mod$model.geographic_extent[[i]]))))
                  z2 <- c(z1[4], z1[1], z1[2], z1[1], z1[2], 
                    z1[3], z1[4], z1[3])
                  x1 <- as.data.frame(matrix(z2, ncol = 2, byrow = TRUE))
                  x1 <- cbind(x1[2], x1[1])
                  p1 <- Polygon(x1)
                  ps1 <- Polygons(list(p1), 1)
                  ecobase_poly2[[i]] <- st_as_sf(SpatialPolygons(list(ps1)))
                }
                ecobase_poly2[[i]] <- ecobase_poly2[[i]]
            }
            for (i in 1:length(ecobase_poly2)) {
                if (!any(class(ecobase_poly2[[i]]) == "sf")) {
                  t2 <- ecobase_poly2[[i]]
                  t3 <- st_cast(t2, to = "POLYGON")
                  ecobase_poly2[[i]] <- st_as_sf(as(st_zm(st_geometry(t3)), 
                    "Spatial"))
                }
                else message("Ok!")
            }
            table2 <- as.data.frame(cbind(1:length(ecobase_poly2), 
                as.numeric(as.character(liste_mod$model.model_number))))
            names(table2) <- c("row", "id")
            id_selected <- as.numeric(as.character(l3$model.model_number))
            rows_selected <- c()
            for (i in 1:length(id_selected)) {
                rows_selected[i] <- as.numeric(table2[table2["id"] == 
                  id_selected[i], ][1])
            }
            spatial1 <- ecobase_poly2[rows_selected]
        }
    }
    if (db == "mg") {
        message("####################### MANGAL DATABASE #######################\n\n")
        message("Fetching datasets from the Mangal website! \nThis operation might take a long time!")
        if (is.null(mangal_types)) 
            mangal_types <- c("predation", "herbivory")
        if ("all" %in% mangal_types) {
            mangal_types <- c("competition", "predation", "herbivory", 
                "amensalism", "neutralism", "commensalism", "mutualism", 
                "parasitism", "symbiosis", "scavenger", "detritivore")
            message("You are downloading the all types of interactions in the mangal database:\n              competition, predation, herbivory, amensalism, neutralism, commensalism,\n              mutualism, parasitism, symbiosis, scavenger, detritivore")
        }
        else mangal_types <- mangal_types
        ntypes <- length(mangal_types)
        net_info <- list()
        type_info <- c()
        for (i in 1:ntypes) {
            message(paste0("\n\nFetching information from interactions of the type ", 
                "'", mangal_types[i], "'!"))
            df_inter <- search_interactions(type = mangal_types[i], 
                verbose = TRUE)
            if (nrow(df_inter) > 0) 
                fwlist1 <- get_collection(df_inter, verbose = TRUE)
            if (nrow(df_inter) > 0) 
                net_info <- c(net_info, fwlist1)
            if (nrow(df_inter) > 0) 
                fwlist2 <- rmangal::as.igraph(fwlist1)
            if (nrow(df_inter) > 0) 
                type_info <- c(type_info, rep(mangal_types[i], 
                  length(fwlist2)))
            if (nrow(df_inter) > 0) 
                fwlist <- c(fwlist, fwlist2)
        }
        for (i in 1:length(fwlist)) {
            fw2 <- fwlist[[i]]
            fw3 <- igraph::as_data_frame(fw2, what = "both")
            id_name <- fw3$vertices[, 1:2]
            for (j in 1:nrow(id_name)) {
                node_name <- (paste0(id_name$original_name[j], 
                  "_", id_name$name[j]))
                if (grepl(":", node_name, fixed = TRUE)) {
                  node_name <- tail(strsplit(node_name, ": "))[[1]]
                  id_name[j, 2] <- node_name[2]
                }
                else id_name[j, 2] <- node_name
            }
            id_edges <- fw3$edges[, 1:3]
            int_matrix <- as.data.frame(matrix(ncol = nrow(id_name), 
                nrow = nrow(id_name)))
            colnames(int_matrix) <- id_name$original_name
            rownames(int_matrix) <- id_name$original_name
            for (a in 1:nrow(id_edges)) {
                edge1 <- as.numeric(id_edges[a, 1:2])
                name1 <- id_name[as.character(edge1[1]), ][, 
                  2]
                name2 <- id_name[as.character(edge1[2]), ][, 
                  2]
                int_matrix[name1, name2] <- 1
            }
            int_matrix[is.na(int_matrix)] <- 0
            fwlist[[i]] <- int_matrix
        }
        references <- as.data.frame(matrix(ncol = 6))
        names(references) <- c("Dataset ID", "Type of interaction", 
            "Original ID", "first_author", "year", "DOI")
        message("Fetching references!")
        for (j in 1:length(net_info)) {
            dataset_id <- net_info[[j]]$dataset$dataset_id
            first_author <- net_info[[j]]$reference$first_author
            year_mng <- as.numeric(net_info[[j]]$reference$year)
            doi_mng <- net_info[[j]]$reference$doi
            references[j, 3] <- dataset_id
            references[j, 4] <- first_author
            references[j, 5] <- year_mng
            references[j, 6] <- doi_mng
            references <- references[order(references$`Dataset ID`), 
                ]
            rownames(references) <- 1:nrow(references)
        }
        references[, 1] <- paste0("mg_", 1:nrow(references))
        references[, 2] <- type_info
        names(fwlist) <- paste0("mg_", references[, 1])
        if (spatial == TRUE) {
            spatial1 <- as.data.frame(matrix(ncol = 4))
            names(spatial1) <- c("Dataset ID", "first_author", 
                "lat", "long")
            message("Fetching coordinates!")
            for (z in 1:length(net_info)) {
                dataset_id <- net_info[[z]]$dataset$dataset_id
                lat_mng <- net_info[[z]]$network$geom_lat
                long_mng <- net_info[[z]]$network$geom_lon
                first_author <- net_info[[z]]$reference$first_author
                if (length(unlist(lat_mng)) > 1) {
                  spatial2 <- as.data.frame(matrix(ncol = 4))
                  names(spatial2) <- c("Dataset ID", "first_author", 
                    "long", "lat")
                  for (b in 1:length(unlist(lat_mng))) {
                    spatial2[b, 3] <- long_mng[[1]][b]
                    spatial2[b, 4] <- lat_mng[[1]][b]
                  }
                  spatial2[, 1] <- dataset_id
                  spatial2[, 2] <- first_author
                  spatial1 <- rbind(spatial1, spatial2)
                }
                spatial1[z, 1] <- dataset_id
                spatial1[z, 2] <- first_author
                if (length(unlist(lat_mng)) == 1) 
                  spatial1[z, 3] <- lat_mng
                if (length(unlist(lat_mng)) == 1) 
                  spatial1[z, 4] <- long_mng
            }
            spatial1 <- spatial1[order(spatial1$`Dataset ID`), 
                ]
            rownames(spatial1) <- 1:nrow(spatial1)
        }
    }
    message(paste0("DONE! \n\nOverall the list stores ", length(fwlist), 
        " datasets!"))
    master_list <- list()
    master_list[["int_matrix"]] <- fwlist
    if (ecosyst == TRUE) {
        master_list[["ecosystem"]] <- ecosystem
        message("\n Additional element in the results: \n\n The vector with information on the ecosystems.")
    }
    if (ref == TRUE) {
        master_list[["references"]] <- references
        message("Additional element in the results! \nA data frame with information on the references.")
    }
    if (spatial == TRUE) {
        master_list[["spatial_info"]] <- spatial1
        message("\n Additional element in the results: \n\n Spatial information was added.")
    }
    if (code == TRUE) {
        if (db == "gw") 
            master_list[["code"]] <- files_gw
        if (db == "wl") 
            master_list[["code"]] <- files_wl
        if (db == "mg") 
            master_list[["code"]] <- references[, 1]
        message("Added food web code information.")
    }
    if (length(master_list) == 1) 
        return(fwlist)
    if (length(master_list) != 1) 
        return(master_list)
    message("####################### DONE! #######################")
}

robustness<-function (fw_to_attack, probs_of_fw, alpha1){
    isolates <- function(g) {
        return(which(degree(g) == 0) - 1)
    }
    fw_nodes <- V(fw_to_attack)$name
    n_species <- length(fw_nodes)
    i_output_list <- list()
    for (j in 1:ncol(probs_of_fw)) {
        output_nodes_and_links <- data.frame(matrix(ncol = n_species, 
            nrow = 5))
        colnames(output_nodes_and_links) <- paste0("del_species_", 
            1:(n_species - 1))
        rownames(output_nodes_and_links) <- c("nodes", "links", 
            "secondary_extinctions", "%_extinctions", "n_primary_extinctions")
        probs_i <- probs_of_fw[, j]
        n_links_original <- gsize(fw_to_attack)
        perc_ext <- 0
        for (z in 1:length(fw_nodes)) {
            node_to_kill <- sample(x = fw_nodes, size = z, replace = FALSE, 
                prob = probs_i)
            fw_resulting_from_attack <- delete_vertices(fw_to_attack, 
                node_to_kill)
            secondary_extinctions <- isolates(fw_resulting_from_attack)
            fw_resulting_from_attack <- igraph::delete_vertices(fw_resulting_from_attack, 
                names(secondary_extinctions))
            if (gsize(fw_resulting_from_attack) == 0 && gorder(fw_resulting_from_attack) == 
                0) {
                output_nodes_and_links[1, z] <- NA
                output_nodes_and_links[2, z] <- NA
                output_nodes_and_links[3, z] <- NA
                output_nodes_and_links[4, z] <- NA
                output_nodes_and_links[5, z] <- NA
            }
            if (gorder(fw_resulting_from_attack) != 0) {
                nodes_in_original_fw <- n_species
                nodes_in_fw_resulting_from_attack <- gorder(fw_resulting_from_attack)
                lost_nodes <- nodes_in_original_fw - nodes_in_fw_resulting_from_attack
                number_of_secondary_extinctions <- length(secondary_extinctions)
                links_in_original_fw <- n_links_original
                links_in_fw_resulting_from_attack <- gsize(fw_resulting_from_attack)
                lost_links <- links_in_original_fw - links_in_fw_resulting_from_attack
                perc_ext <- ((n_species - nodes_in_fw_resulting_from_attack) * 
                  100)/n_species
                output_nodes_and_links[1, z] <- nodes_in_fw_resulting_from_attack
                output_nodes_and_links[2, z] <- links_in_fw_resulting_from_attack
                output_nodes_and_links[3, z] <- number_of_secondary_extinctions
                output_nodes_and_links[4, z] <- perc_ext
                output_nodes_and_links[5, z] <- z
            }
        }
        i_output_list[[j]] <- output_nodes_and_links
    }
    R_alpha1 <- NA
    if (alpha1 == 50) 
        prop_species <- n_species/2
    if (alpha1 == 100) 
        prop_species <- n_species
    for (x in 1:length(i_output_list)) {
        df1 <- i_output_list[[x]]
        df1[4, ] <- round(as.numeric(df1[4, ]))
        if (any(as.numeric(df1[4, ]) <= alpha1, na.rm = TRUE)) {
            col1 <- max(which(df1[4, ] <= alpha1))
            R_alpha1[x] <- df1[5, col1]/n_species
        }
        else {
            R_alpha1[x] <- 0
        }
    }
    return(list(Simulation_results = i_output_list, Ralpha = R_alpha1))
}
