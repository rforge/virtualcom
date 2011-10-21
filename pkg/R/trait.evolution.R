trait.evolution <- function(branchingRate=0.1, Nleaves=100, Ninv=1, minOpt=0, maxOpt=100, which.evolution.model="BM", mytheta=1, mysigma=.01){
	# branchingRate <- 0.1	
	# Nleaves <- 30	# nb of species
	# Ninv <- 10 # nb of invasive species
	# minOpt <- 0 # minimal optimal niche value
	# maxOpt <- 100 # maximal optimal niche value
	# mytheta <- 1 # strength of attractor in OU process
	# to test: branchingRate=0.1; Nleaves=30; Ninv=1; minOpt=0; maxOpt=100; which.evolution.model="BM"; mytheta=1; mysigma=.01
	
	# Brownian motion of tree
	myTree <- treedata(birthdeath.tree(branchingRate, 0, taxa.stop=(Nleaves + 1)), data.frame(rnorm(Nleaves)), warnings=FALSE)$phy
	# plot(myTree)
	# change naming of tips so that tip labels and tip numbers in edge are the same and are ordered from 1 to Nleaves
	myTree$edge[myTree$edge <= length(myTree$tip.label)] <- 1: length(myTree$tip.label)
	myTree$tip.label <- 1: length(myTree$tip.label)

 	# different models of trait evolution
 	logist <- function(x, lower=0, upper=100){((x-min(x))*(upper - lower)/abs(max(x)-min(x))) + lower}  # To rescale the trait value between selected boundaries
	if (which.evolution.model == "BM") myTrait <- rTraitCont(myTree, model="BM", sigma=mysigma,  root.value=50)
	if (which.evolution.model == "OU"){
        nGroups = round(Nleaves/10)
        thetaLeaves <- seq(0, 100, length.out=nGroups)
        attractors <- get.attractors(phy=myTree, theta=thetaLeaves, ngroups=nGroups, method="single", plotting=FALSE, randomize=TRUE)
        attractors$attractor.r <- rnorm(nrow(attractors), attractors$attractor, 1)
        if (nGroups == 1) { mymu <- 0} else {
           mymu <- attractors[match(myTree$edge[,2], attractors$node.index), "attractor.r"]
        }
        myTrait <- mymu[myTree$edge[,2] < (1+length(myTree$tip.label))]
        names(myTrait) <- myTree$edge[,2][myTree$edge[,2] < (1+length(myTree$tip.label))]        
        # myTrait <- rTraitCont(myTree,model="OU",theta=mymu,alpha=0.1,sigma=mysigma) # TODO: is not working at the moment, why?
        # phylosignal(myTrait,myTree);summary(myTrait)
        # myTrait <- evolve.trait(myTree, trait.mode="OU", x.root=50, sigma=10, mu=mymu, theta=mytheta, show.plots.trait = FALSE, use.color=TRUE, main=paste("theta=",mytheta,sep=""))[[1]] # this works but is too slow
  	}
  	if (which.evolution.model == "delta"){
        myTree_t <- deltaTree(myTree, 0.01)
        myTrait <- rTraitCont(myTree_t, model="BM", sigma=mysigma,  root.value=50)
  	}
  	# rescale traits to make them comparable among different trait evolution models
  	# myTraitold <- myTrait
  	myTrait <- logist(myTrait, lower=minOpt, upper=maxOpt)
	# plot(myTraitold, myTrait)
 	opt <- seq(minOpt,maxOpt,length.out=Nleaves) # optimum value for all species (equal spacing between the species)
 	attra <- as.data.frame(cbind(SpeciesID = myTree$tip.label, niche_evol=myTrait, niche_equi=opt[as.vector(rank(myTrait))]))  
 	rownames(attra) <- attra$SpeciesID
  	# test: sum(rank(attra$niche_evol) != rank(attra$niche_equi))
  	# PS of evolved and equidistant trait: phylosignal(attra$niche_evol, myTree);  phylosignal(attra$niche_equi, myTree) 

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
		attra[names(sort(cophenetic.phylo(myTree)[,uu])[1:10]),"UnderInv"] <- 1   
	}
	attra <- attra[order(attra$SpeciesID),]       
	return(list(attra, myTree))
}
   
