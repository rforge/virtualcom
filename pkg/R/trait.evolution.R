trait.evolution <- function(branchingRate = 0.1, Nleaves = 100, Ninv = 1, minOpt = 0, maxOpt = 100, which.evolution.model = "BM", rescale=TRUE, mysigma = 0.01, tree.shape.param = NA, OUwie.sigma.sq=NA, OUwie.theta=NA, OUwie.alpha=NA) {
    # Brownian motion of tree
    myTree <- sim.bdtree(branchingRate, 0, stop = "taxa", n = Nleaves)
    # change naming of tips so that tip labels and tip numbers in edge are the same and are ordered from 1 to Nleaves
    myTree$edge[myTree$edge <= length(myTree$tip.label)] <- 1:length(myTree$tip.label)
    myTree$tip.label <- 1:length(myTree$tip.label)
    
    # different models of trait evolution
    if (which.evolution.model == "BM") {
        myTrait <- rTraitCont(myTree, model = which.evolution.model, sigma = mysigma, root.value = 0)
    }
    if (which.evolution.model == "deltaTree") {
        myTree_tmp <- rescale(myTree, model = "delta", tree.shape.param)
        # to capture the cases where rescale transformes branch lengths into NA or zero or negative values
        if (any(is.na(myTree_tmp$edge.length))) {
            which <- which(is.na(myTree_tmp$edge.length))
            myTree_tmp$edge.length[which] <- min(myTree_tmp$edge.length, na.rm = T)
            myTree_tmp$edge.length[-which] <- myTree_tmp$edge.length[-which] + myTree_tmp$edge.length[which]
        }
        myTrait <- rTraitCont(myTree_tmp, model = "BM", sigma = mysigma, root.value = 0)
        names(myTrait) <- 1:length(myTrait)
    }
    if (which.evolution.model == "lambdaTree") {
        myTree_tmp <- rescale(myTree, model = "lambda", tree.shape.param)
        myTrait <- rTraitCont(myTree_tmp, model = "BM", sigma = mysigma, root.value = 0)
        names(myTrait) <- 1:length(myTrait)
    }
    if (which.evolution.model == "kappaTree") {
        myTree_tmp <- rescale(myTree, model = "kappa", tree.shape.param)
        myTrait <- rTraitCont(myTree_tmp, model = "BM", sigma = mysigma, root.value = 0)
        names(myTrait) <- 1:length(myTrait)
    }
    if (which.evolution.model == "OUwie.sim") {
        myTree_tmp <- myTree
        datas <- get.attractors(phy=myTree_tmp, theta=OUwie.theta, ngroups=length(OUwie.theta), method="single", plotting=FALSE, randomize=TRUE)
        myTree_tmp$node.label <- datas$attractor[(length(myTree_tmp$tip.label)+1) : nrow(datas)]
        if(length(levels(factor(myTree_tmp$node.label))) != length(OUwie.theta)) warning("OUwie.sim does not work for single-tip attractors")
         myTrait <- OUwie.sim(phy=myTree_tmp, data=datas[1:length(myTree_tmp$tip.label),], sigma.sq = OUwie.sigma.sq, theta0 = 0, theta=OUwie.theta, alpha=OUwie.alpha)[,3]
        names(myTrait) <- 1:length(myTrait)
    }
    
    # To rescale the trait value between selected boundaries
    logist <- function(x, lower = 0, upper = 100) {
        ((x - min(x)) * (upper - lower)/abs(max(x) - min(x))) + lower
    }
    if(rescale==TRUE) myTrait <- logist(myTrait, lower = minOpt, upper = maxOpt)
    opt <- seq(minOpt, maxOpt, length.out = Nleaves)  # optimum value for all species (equal spacing between the species)
    attra <- as.data.frame(cbind(SpeciesID = myTree$tip.label, niche_evol = myTrait))
    rownames(attra) <- attra$SpeciesID
    
    # Identifying the invaders in the tree:
    if (Ninv == 1) {
        attra$RandInv <- 0
        attra$RandInv[sample(1:nrow(attra), 1)] <- 1
    }
    if (Ninv > 1) {
        # Invaders are overdispersed:
        myDist <- as.dist(as.data.frame(cophenetic.phylo(myTree)), upper = TRUE, diag = TRUE)  # calculate distance matrix for tree
        myCluster <- hclust(myDist, method = "single")  # plot(myCluster) ; make the clusters
        myGroups <- cutree(myCluster, k = Ninv)  # cut the tree
        attra$Group <- myGroups
        attra$OverInv <- 0
        inva <- aggregate(attra$SpeciesID, list(attra$Group), function(x) {
            ii <- sample(as.character(x), 1)
            ii
        })
        attra[as.character(inva[, 2]), "OverInv"] <- 1
        # Invaders are random-dispersed:
        attra$RandInv <- sample(attra$OverInv, length(attra$OverInv))
        # Invaders are underdispersed:
        attra$UnderInv <- 0
        uu <- as.character(sample(as.character(attra$SpeciesID), 1))
        attra[names(sort(cophenetic.phylo(myTree)[, uu])[1:Ninv]), "UnderInv"] <- 1
    }
    attra <- attra[order(attra$SpeciesID), ]
    return(list(attra, myTree))
} 