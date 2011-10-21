param.test.inv <- function(parameters, ...){
  	# to test: parameters=as.numeric(param[1,]); names(parameters) <- names(param[1,]);
  	output <- list()
  	
    #-------------------------------------
    # get input parameters
    #-------------------------------------  
    my.names <- names(parameters)
  	parameters <- as.numeric(parameters); 
    names(parameters) <- my.names
    # general
  	mytimesteps <- parameters["timesteps"]		# timesteps in each simulation to get equilibrium conditions
  	invTime <- parameters["invasionTime"]		# how many years can invaders try to invade
 	myNulRep <- parameters["NulRep"]			# how many repetitions in null models
	# species pool
  	myN <- parameters["N"] 						# size of species pool
  	myNinv <- parameters["Ninv"] 				# number of invaders in species pool
  	myniche.breath <- parameters["niche.breath"]# widths of species niches
 	evolModel <- ifelse(parameters["evolModel"]==1, "BM", ifelse(parameters["evolModel"]==2, "OU", "delta")) # choice of evolutionary model which determines phylogenetic signal
  	myNiches <- parameters["NicheVal"]			# take the evolved traits or equidistant traits 
  	myphySig <- parameters["phySig"] 			# min value of phylogenetic signal in species pool
  	# communities
  	mynbCom <- parameters["nbCom"]				# number of sinulated communities
  	mycapa <- parameters["capacity"] 			# carrying capacity of communities
 	myenvOptimum <- parameters["envOptimum"]	# optimal value of environmental conditions for communities
  	mycstE <- parameters["cstE"]				# strength of environmental filter
  	mycstB <- parameters["cstB"]				# strength of biotic filter
  	mycstSR <- parameters["cstSR"]				# strength of local seed reproduction
 	mycomp <- ifelse(parameters["comp"]==1, "NO", ifelse(parameters["comp"]==2,"NN", ifelse(parameters["comp"]==3,"NO_cut", NA))) # different descriptions of biotic filter (competition)
	myNN <- parameters["NN"]					# how many neighbours included for biotic filter (comeptition)
  	output$parameter <- parameters
  
    #-------------------------------------
    # Species pool: phylogenetic tree, trait values and invaders in the tree
    #-------------------------------------  
   	# Let trees and traits evolve
   	species <- trait.evolution(branchingRate=0.1, Nleaves=myN, Ninv=myNinv, which.evolution.model=evolModel, mytheta=1) 
   	if(myNiches == 1) myopt <- species[[1]]$niche_evol else myopt <- species[[1]]$niche_equi
  	names(myopt) <- species[[1]]$SpeciesID  
  	
  	# Redo if phylogenetic signal should not be below a certain value:
  	if(!is.na(myphySig)){ 
  		## If K < myphySig recreate the tree and traits again until K >= myphySig:
  		while(phylosignal(myopt, species[[2]])$K < myphySig){   
  			species <- trait.evolution(branchingRate=0.1, Nleaves=myN, Ninv=myNinv, which.evolution.model=evolModel, mytheta=1) 
  			if(myNiches == 1) myopt <- species[[1]]$niche_evol else myopt <- species[[1]]$niche_equi 
  			names(myopt) <- species[[1]]$SpeciesID  
  		}
  	}
  	names(species) <-  c("traits", "phylo")
    species$K <- phylosignal(myopt, species[[2]])$K
    species$invader <- list(ID=species[[1]][species[[1]]$Inv==1,"SpeciesID"], opt=myopt[species[[1]][species[[1]]$Inv==1,"SpeciesID"]], distToOpt= abs(myenvOptimum - myopt[species[[1]][species[[1]]$Inv==1,"SpeciesID"]]))   
    output$pool <- species
    
    #-------------------------------------
    # Native community assembly 
    #-------------------------------------
    # Get the optima for the native species only 
    myopt_Nat <- myopt[species[[1]]$Inv==0]	
  	names(myopt_Nat) <- species[[1]][species[[1]]$Inv==0,"SpeciesID"]  

    # Create the list for storing results for the native community
  	l.com.nat <- list()
  	# Initialization: creates the siteXspecies matrix
  	com.matrx <- matrix(0,nrow=mynbCom,ncol=myN,dimnames=list(c(1:mynbCom),c(species[[1]]$SpeciesID)))
  	
    # Community assembly:random
  	if(mycstB==0 && mycstE==0 && mycstSR==0){
  		com.matrx2 <- com.matrx
  		for (i in 1:mynbCom){
  			outspcC <- as.character(sample(species[[1]][species[[1]]$Inv==0,"SpeciesID"], size=mycapa,replace=TRUE)) # get IDs of species
  	 		tt <- table(outspcC) 
       		vv <- as.vector(tt)  
       		names(vv) <- names(tt) # vector of present species
       		l.com.nat$com[[i]] <- list(communities=outspcC, spcXs=t(as.matrix(vv))) 
          	com.matrx2[i,] <- spsXsites(i=i, Nsp=myN, my.matrix=com.matrx2, my.list=l.com.nat$com) # species times site matrix for present and absent species                    
     	}
    } else {
   	# Community assembly: with abiotic filters and biotic interactions and reproduction/seed rain
    	com.matrx2 <- com.matrx
    	for (i in 1:mynbCom){
   			l.com.nat$com[[i]] <- tamaure(niche.breath=myniche.breath, opt=myopt_Nat, env=myenvOptimum, cstE=mycstE, cstB=mycstB, cstSR=mycstSR, years=mytimesteps, kcapa=mycapa, communityIn=NA, plot=FALSE, mypar=c(2,5), comp=mycomp, NN=myNN)
          	com.matrx2[i,] <- spsXsites(i=i, Nsp=myN, my.matrix=com.matrx2, my.list=l.com.nat$com) # species times site matrix for present and absent species 
  		}
    }
    output$native$communities <- l.com.nat$Site_X_Species <- com.matrx2  
      
    # Get indices and null models for diversity:
 	Cophi <- cophenetic(species[[2]])
 	Cofun <- as.matrix(dist(myopt,diag=TRUE,upper=TRUE))
    output$native$indices <- div.param.native(spSite=output$native$communities, phy=Cophi, fun=Cofun, nrep=100, null.model = "taxa.labels") # zNULL = NaN when sdNULL=0				  
    #-------------------------------------
    # Invaded community assembly
    #-------------------------------------
    # Create the list for storing results for the invaded community 
    l.com.inv <- list()
    # list of native community used to start invasion process
    finalCom <- lapply(l.com.nat$com, function(x) x$communities) 
     
    # Invasion assembly: random 
  	if(mycstB==0 && mycstE==0 && mycstSR==0){
		com.matrx3 <- com.matrx
  		for (i in 1:mynbCom){
  			outspcC <- as.character(sample(species[[1]]$SpeciesID, size=mycapa,replace=TRUE))  # get IDs of species
  	 		tt <- table(outspcC)
  	 		vv <- as.vector(tt) 
       		names(vv) <- names(tt)  # vector of present species
       		l.com.inv$com[[i]] <- list(communities=outspcC, spcXs=t(as.matrix(vv))) 
          	com.matrx3[i,] <- spsXsites(i=i, Nsp=myN, my.matrix=com.matrx3, my.list=l.com.inv$com) # species times site matrix for present and absent species  
      	}
      	output$invasion$communities <- l.com.inv$Site_X_Species <- com.matrx3
  	  	output$invasion$process <- NA      
  	} else {
  	# Invasion assembly: with abiotic filters and biotic interactions and reproduction/seed rain
  		com.matrx3 <- com.matrx
  		for (i in 1:mynbCom){
  	  		invaded  <- vector()			# to count the invasion time
  	  		com.inv <- list()
  	   		lastYear <- finalCom[[i]]
          	for(ti in 1:invTime){			# to record the invasion dynamics
          		com.inv[[ti]] <- tamaure(niche.breath=myniche.breath, opt=myopt, env=myenvOptimum, cstE=mycstE, cstB=mycstB, cstSR=mycstSR, years=1, kcapa=mycapa, communityIn=lastYear, plot=FALSE, mypar=c(2,5), comp=mycomp, NN=myNN)
            	lastYear <- com.inv[[ti]]$communities
      			invaded[ti] <- ifelse(as.character(species$invader$ID) %in% names(table(lastYear)), table(lastYear)[as.character(species$invader$ID)], 0)
  			}
  			l.com.inv$com[[i]] <- com.inv[[invTime]]
  			l.com.inv$process[[i]] <- list(inv=invaded, com=com.inv)
          	com.matrx3[i,] <- spsXsites(i=i, Nsp=myN, my.matrix=com.matrx3, my.list=l.com.inv$com) # species times site matrix for present and absent species
  	  	}
  	  	output$invasion$communities <- l.com.inv$Site_X_Species <- com.matrx3
  	  	output$invasion$process <- t(sapply(l.com.inv$process, function(x) x$inv))	# the community years are in line and the years are in colomne
    }		
   
    # Get indices for similarity between inv and nat species and null models:
 	Cophi <- cophenetic(species[[2]])
 	Cofun <- as.matrix(dist(myopt,diag=TRUE,upper=TRUE))
 	mysamp <- as.data.frame(output$invasion$communities)
  	colnames(mysamp) <- paste("sp",colnames(mysamp),sep="")
  	colnames(Cophi) <- rownames(Cophi) <- paste("sp",colnames(Cophi),sep="")
  	colnames(Cofun) <- rownames(Cofun) <- paste("sp",colnames(Cofun),sep="")
    # make sure order is correct 	  
    Cophi <- Cophi[colnames(mysamp), colnames(mysamp)]
    Cofun <- Cofun[colnames(mysamp), colnames(mysamp)]
    myinva <- paste("sp", species$invader$ID,sep="")
    
    if (is.na(output$invasion$process) && sum(output$invasion$communities[,species$invader$ID])==0){
    	output$invasion$indices <- NA} else {
    	if (is.na(output$invasion$process) && sum(output$invasion$communities[,species$invader$ID])!=0) output$invasion$indices <- div.param.invasion(spSite=mysamp, phy=Cophi, fun=Cofun, nrep=myNulRep, inva=myinva, null.model = "native_inv")
    	if(!is.na(output$invasion$process) && sum(output$invasion$process[,invTime]) != 0) output$invasion$indices <- div.param.invasion(spSite=mysamp, phy=Cophi, fun=Cofun, nrep=myNulRep, inva=myinva, null.model = "native_inv")
    	if(!is.na(output$invasion$process) && sum(output$invasion$process[,invTime]) == 0) output$invasion$indices <- NA
    }	
  	return(output)
}   