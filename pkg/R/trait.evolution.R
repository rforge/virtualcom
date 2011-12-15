trait.evolution <- function(branchingRate=0.1, Nleaves=100, Ninv=1, minOpt=0, maxOpt=100, which.evolution.model="BM", mysigma=.01, extraTreeParam=NA){
	# to test: branchingRate=0.1; Nleaves=30; Ninv=1; minOpt=0; maxOpt=100; which.evolution.model="BM"; mysigma=.01; extraTreeParam=NA
  # to test: branchingRate=0.1; Nleaves=n.species.pool; Ninv=n.invader.pool; which.evolution.model=evol.model; extraTreeParam=evol.model.param
  
	# Brownian motion of tree
	# myTree <- treedata(birthdeath.tree(branchingRate, 0, taxa.stop=(Nleaves + 1)), data.frame(rnorm(Nleaves)), warnings=FALSE)$phy
	myTree <- drop.tip(birthdeath.tree(b=branchingRate, d=0, taxa.stop=(Nleaves + 1)), as.character((Nleaves + 1)))
  
  # change naming of tips so that tip labels and tip numbers in edge are the same and are ordered from 1 to Nleaves
	myTree$edge[myTree$edge <= length(myTree$tip.label)] <- 1: length(myTree$tip.label)
	myTree$tip.label <- 1: length(myTree$tip.label)

 	# different models of trait evolution  
   if (which.evolution.model == "BM"){
       myTrait <- rTraitCont(myTree, model=which.evolution.model, sigma=mysigma,  root.value=0)
       if (!is.na(extraTreeParam)){
          traitA <- sample(myTrait)
          names(traitA) <- 1:length(traitA)
          myTrait <- sqrt(extraTreeParam) * myTrait + sqrt(1-extraTreeParam) * traitA
      }
   }
   if (which.evolution.model == "deltaTree"){
       myTree_tmp <- deltaTree(phy=myTree, delta=extraTreeParam, rescale = T)
       myTrait <- rTraitCont(myTree_tmp, model="BM", sigma=mysigma,  root.value=0)
       names(myTrait) <- 1:length(myTrait)
   }
   if (which.evolution.model == "lambdaTree"){
       myTree_tmp <- lambdaTree(phy=myTree, lambda=extraTreeParam)
       myTrait <- rTraitCont(myTree_tmp, model="BM", sigma=mysigma,  root.value=0)
       names(myTrait) <- 1:length(myTrait)
   }
   if (which.evolution.model == "kappaTree"){
       myTree_tmp <- kappaTree(phy=myTree, kappa=extraTreeParam)
       myTrait <- rTraitCont(myTree_tmp, model="BM", sigma=mysigma,  root.value=0)
       names(myTrait) <- 1:length(myTrait)
   }
   if (which.evolution.model == "ouTree"){
       myTree_tmp <- ouTree(phy=myTree, alpha=extraTreeParam)
       myTrait <- rTraitCont(myTree_tmp, model="BM", sigma=mysigma,  root.value=0)
       names(myTrait) <- 1:length(myTrait)
   }   
 	
 	 # To rescale the trait value between selected boundaries
 	 logist <- function(x, lower=0, upper=100){((x-min(x))*(upper - lower)/abs(max(x)-min(x))) + lower}  

   # rescale traits to make them comparable among different trait evolution models
 	 myTrait <- logist(myTrait, lower=minOpt, upper=maxOpt)
 	 opt <- seq(minOpt, maxOpt, length.out=Nleaves) # optimum value for all species (equal spacing between the species)
 	 attra <- as.data.frame(cbind(SpeciesID = myTree$tip.label, niche_evol=myTrait, niche_equi=opt[as.vector(rank(myTrait))]))  
 	 rownames(attra) <- attra$SpeciesID

	 #  Identifying the invaders in the tree:
	 if(Ninv==1){
      attra$Inv <- 0
      attra$Inv[sample(1:nrow(attra),1)] <- 1	
	 }
	 if(Ninv > 1){
	  	# Invaders are overdispersed:
		  myDist <- as.dist(as.data.frame(cophenetic.phylo(myTree)),upper=TRUE,diag=TRUE) # calculate distance matrix for tree
		  myCluster <- hclust(myDist, method="single") # plot(myCluster) ; make the clusters
		  myGroups <- cutree(myCluster, k=Ninv)	# cut the tree
		  attra$Group <- myGroups
		  attra$OverInv <- 0
		  inva <- aggregate(attra$SpeciesID, list(attra$Group), function(x) {ii <- sample(as.character(x),1); ii})
		  attra[as.character(inva[,2]),"OverInv"] <- 1

    	# Invaders are random-dispersed:
		  attra$RandInv <- sample(attra$OverInv,length(attra$OverInv))

    	# Invaders are underdispersed:
		  attra$UnderInv <- 0
		  uu <- as.character(sample(as.character(attra$SpeciesID),1))      
		  attra[names(sort(cophenetic.phylo(myTree)[,uu])[1:Ninv]),"UnderInv"] <- 1   
	}
	attra <- attra[order(attra$SpeciesID),]       
	return(list(attra, myTree))
}
   
